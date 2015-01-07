#!/usr/bin/env python

from __future__ import division
import numpy as np
import squmfit
from squmfit import Fit, Model
from matplotlib import pyplot as pl
from collections import namedtuple

class Aniso(object):
    """ Anisotropy detection channels """
    def __init__(self, par, perp):
        self.par = par
        self.perp = perp

    def map(self, f):
        return Aniso(f(self.par), f(self.perp))

def exponential(t, rate, amplitude):
    """ Note that this should contain a prefactor of `rate` that has
    been omitted to minimize parameter covariance """
    return amplitude * np.exp(-t * rate)

ExponentialModel = Model(exponential)

class ConvolvedModel(squmfit.Expr):
    """ A model convolved with a response function """
    def __init__(self, response, period, model, offset=0):
        """
        Create a convolved model

        :param response: The instrument response
        :param period: The period of the instrument response
        :param model: The model we are trying to evaluate
        :param offset: The temporal offset between the response and
        the measured curve (in bins)
        """
        assert isinstance(model, squmfit.Expr)
        self.response = response
        self.period = period
        self.model = model
        self.offset = offset

    def evaluate(self, params, **user_args):
        from scipy.signal import fftconvolve
        from scipy.interpolate import interp1d

        offset = self.offset
        if isinstance(self.offset, squmfit.Expr):
            offset = offset.evaluate(params, **user_args)

        period = self.period
        if isinstance(self.period, squmfit.Expr):
            period = period.evaluate(params, **user_args)

        model = self.model.evaluate(params, **user_args)

        shift = len(model) % period
        periodic_irf_f = interp1d(np.arange(len(self.response)), self.response,
                                  assume_sorted=True)
        ts = np.arange(10 * len(self.response)) + shift + offset
        ts %= period
        interpolated_irf = periodic_irf_f(ts)

        def pad(arr, n):
            m = arr.shape[0]
            if m < n:
                return np.hstack([arr, np.zeros(n-m)])
            else:
                return arr

        #model = pad(model, n)
        a = fftconvolve(interpolated_irf, model, 'same')
        n_periods = 10
        n = round(n_periods * period)
        return a[n:n + len(model)]

    def parameters(self):
        return self.model.parameters()

def estimate_rep_rate(irf):
    """ Estimate excitation repetition rate (measured in bins) from IRF """
    middle = (np.max(irf) - np.min(irf)) / 2 + np.min(irf)
    idxs, = np.nonzero(np.logical_and(irf[:-1] < middle, middle < irf[1:]))
    a,b = sorted(idxs)[:2]
    period = b - a
    debug = False
    if debug:
        pl.plot(irf)
        pl.axhline(middle, c='k')
        pl.axvline(a, c='g')
        pl.axvline(b, c='g')
        pl.yscale('log')
        pl.show()
    return period

def normalize_irfs(irfs):
    """
    Subtract background from and normalize IRF

    :type irfs: Aniso of histograms
    :rtype: Aniso of histograms
    """
    def background_subtract(irf):
        bg = np.median(irf)
        print('IRF background = %1.2f' % bg)
        return irf - bg

    irfs = irfs.map(background_subtract)

    # Fix normalization of IRF
    # We need to take care to preserve the relative magnitude of the
    # perpendicular/parallel IRFs
    norm = sum(irfs.par)
    return irfs.map(lambda x: x / norm)

def fit(irfs, corrs, jiffy_ps, exc_period, n_components, periods=1):
    """
    :type irfs: list of `Aniso`s
    :param irfs: IRF histograms
    :type corrs: list of `Aniso`s
    :param corrs: Fluorescence histograms
    :type jiffy_ps: int
    :param jiffy_ps: Bin width in picoseconds
    :type exc_period: int
    :param exc_period: Excitation period measured in bins
    :type n_components: int
    :param n_components: Number of exponential decay components to fit against
    :type periods: int
    :param periods: Number of periods to fit against
    """
    jiffy = jiffy_ps * 1e-12
    n = periods * exc_period
    irfs = irfs.map(lambda x: x[:n])

    # Run the fit
    res = analyze(irfs, corrs, exc_period, n_components, jiffy_ps)
    res = analyze(irfs, corrs, exc_period, n_components, jiffy_ps,
                  free_period=True, params0=res.params)
    return res

def analyze(irfs, corrs, exc_period, n_components, jiffy_ps,
            params0=None, free_period=False):
    fit = Fit()

    offset_par = fit.param('offset-par', 0)
    offset_perp = fit.param('offset-perp', 0)
    period = fit.param('period', exc_period)

    # Build decay model
    rates = []
    for i in range(n_components):
        tau = 1000 + 1000*i
        rate = fit.param('lambda%d' % i, initial=1/tau)
        rates.append(rate)

    # Parameters for anisotropy model
    r0 = fit.param('r0', initial=0.4)
    rate_rot = fit.param('lambda_rot', initial=1/1000)
    imbalance = fit.param('g', initial=1)

    for pair_idx,corr in enumerate(corrs):
        def read_histogram(path):
            # read histogram
            corr = np.genfromtxt(path)[:,1]
            n = len(irfs.par) # FIXME?
            assert len(corr) >= n
            corr = corr[:n] # FIXME?
            return corr

        # generate fluorescence decay model
        par = read_histogram(corr.par)
        perp = read_histogram(corr.perp)
        decay_models = []
        initial_amp = np.max(par) / np.sum(par)
        for comp_idx, rate in enumerate(rates):
            amp = fit.param('c%d_amplitude%d' % (pair_idx, comp_idx), initial=initial_amp)
            decay_models.append(ExponentialModel(rate=rate, amplitude=amp))
        decay_model = sum(decay_models)

        def add_curve(corr, name, rot_model, offset, norm, irf):
            times = jiffy_ps * np.arange(corr.shape[0])

            # generate weights
            weights = np.zeros_like(corr)
            weights[corr != 0] = 1 / np.sqrt(corr[corr != 0])

            # generate model
            convolved = ConvolvedModel(irf,
                                       period if free_period else exc_period,
                                       decay_model * rot_model(times),
                                       offset=offset)
            model = norm * np.sum(corr) * convolved
            fit.add_curve(name, model, corr, weights=weights, t=times)

        add_curve(corr = par,
                  rot_model = lambda times: 1 + 2 * r0 * np.exp(-rate_rot * times),
                  norm = 1,
                  offset = offset_par,
                  name = corr.par+'_par',
                  irf = irfs.par)
        add_curve(corr = perp,
                  rot_model = lambda times: 1 - r0 * np.exp(-rate_rot * times),
                  norm = imbalance,
                  offset = offset_perp,
                  name = corr.perp+'_perp',
                  irf = irfs.perp)

    return fit.fit(params0)

def print_params(p):
    print '  irf period', p['period']
    print '  irf offset (parallel)', p['offset-par']
    print '  irf offset (perpendicular)', p['offset-perp']
    print '  g', p['g']
    print '  r0', p['r0']
    print '  tau_rot', 1/p['lambda_rot']

    for comp_idx in range(args.components):
        rate = p['lambda%d' % comp_idx]
        print '  Component %d' % comp_idx
        print '    tau', 1/rate

    for pair_idx,pair in enumerate(corrs):
        print '  Curve %d' % pair_idx
        for comp_idx in range(args.components):
            rate = p['lambda%d' % comp_idx]
            amp = p['c%d_amplitude%d' % (pair_idx, comp_idx)] / rate
            print '    amplitude%d' % comp_idx, amp

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('corr', metavar='FILE', nargs='+',
                        help="correlation function")
    parser.add_argument('--irf', '-i', metavar='FILE', action='append',
                        help='instrument response function')
    parser.add_argument('--components', '-c', type=int, default=1,
                        help='number of fit components')
    parser.add_argument('--rep-rate', '-r', type=float,
                        help='pulse repetition rate (Hertz)')
    parser.add_argument('--periods', '-p', type=int, default=1,
                        help='how many pulse periods we should fit to')
    parser.add_argument('--output', '-o', type=str,
                        help='where to send output')
    parser.add_argument('--no-offset', action='store_true',
                        help='do not fit temporal offset between data and IRF')
    parser.add_argument('-j', '--jiffy', type=float,
                        help='Bin width')
    args = parser.parse_args()

    corrs = [Aniso(a,b) for a,b in zip(args.corr[::2], args.corr[1::2])]

    irfs = [np.genfromtxt(irf, dtype=None, names='time,counts') for irf in args.irf]
    if len(irfs) != 2:
        raise RuntimeError('Expected two IRFs')
    times = irfs[0]['time']
    irfs = [irf['counts'] for irf in irfs]
    irfs = normalize_irfs(Aniso(irfs[0], irfs[1]))

    # Determine the channel width (jiffy)
    if args.jiffy is not None:
        jiffy_ps = args.jiffy / 1e-12
    else:
        jiffy_ps = (times[1] - times[0]) # in picoseconds

    # Determine the pulse repetition rate
    if args.rep_rate is None:
        period = estimate_rep_rate(irfs.par)
    else:
        period = int(1 / args.rep_rate / jiffy) # period in ticks

    print 'Period', period, 'bins'
    print 'Channel width', jiffy_ps, 'ps'

    res = fit(irfs, corrs, jiffy_ps, period, args.components)

    # Present results
    print
    print 'Initial parameters'
    print_params(res.initial.params)

    print
    print 'Fitted parameters'
    print_params(res.params)

    # Fix covariance
    for comp_idx1 in range(args.components):
        for pair_idx1,_ in enumerate(args.corr):
            for comp_idx2 in range(args.components):
                for pair_idx2,_ in enumerate(args.corr):
                    p1 = 'c%d_amplitude%d' % (pair_idx1, comp_idx1)
                    p2 = 'c%d_amplitude%d' % (pair_idx2, comp_idx2)
                    rate1 = res.params['lambda%d' % comp_idx1]
                    rate2 = res.params['lambda%d' % comp_idx2]
                    #res.covar[p1][p2] *= rate1 * rate2

    print
    print 'Reduced chi-squared'
    for name, curve in res.curves.items():
        print '  %-15s     %1.3g' % (name, curve.reduced_chi_sqr)

    print
    print 'Standard error'
    if res.stderr is not None:
        for param, err in res.stderr.items():
            print '  %-15s     %1.2g' % (param, err)
    else:
        print "  Failed to compute due to flat axis"

    print
    print 'Correlations (coefficients less than 0.2 omitted)'
    if res.correl is not None:
        correls = {(param1,param2): res.correl[param1][param2]
                   for param1 in res.params.keys()
                   for param2 in res.params.keys()
                   if param1 < param2}
        for (p1,p2), c in sorted(correls.items(), key=lambda ((a,b),c): c, reverse=True):
            if abs(c) > 0.2:
                print '  %-15s / %-15s       %1.2f' % (p1, p2, c)
    else:
        print "  Failed to compute due to flat axis"

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    plots = pl.subplot()
    divider = make_axes_locatable(plots)
    residuals = divider.append_axes('bottom', size='10%', pad=0.05)
    color_cycle = pl.rcParams['axes.color_cycle']
    for pair_idx, aniso in enumerate(corrs):
        for name, ch in [(aniso.par, 'par'), (aniso.perp, 'perp')]:
            cres = res.curves['%s_%s' % (name, ch)]
            sym = '+' if ch == 'par' else 'x'
            color = color_cycle[pair_idx % len(color_cycle)]
            times = jiffy * np.arange(cres.fit.shape[0])
            #plots.plot(times, res.initial.curves[name].fit, '+', label='Initial')
            plots.plot(times, cres.curve.data, sym, label='Observed', color=color)
            plots.plot(times, cres.fit, label='Fit', color=color)
            residuals.plot(times, cres.residuals, sym, color=color)

    plots.set_ylabel('number of occurences')
    residuals.set_ylabel('residual')
    residuals.set_xlabel('time (ps)')
    residuals.axhline(0, color='k')
    plots.set_yscale('log')
    plots.legend()
    if args.output is not None:
        pl.savefig(args.output, figsize=(5,5))
    else:
        pl.show()

if __name__ == '__main__':
    main()
