#!/usr/bin/env python

from __future__ import division
import numpy as np
import squmfit
from squmfit import Fit, Model
from matplotlib import pyplot as pl

def exponential(t, rate, amplitude):
    """ Note that this should contain a prefactor of `rate` that has
    been omitted to minimize parameter covariance """
    return amplitude * np.exp(-t * rate)

ExponentialModel = Model(exponential)

class ConvolvedModel(squmfit.Expr):
    def __init__(self, response, period, model, offset=0):
        assert isinstance(model, squmfit.Expr)
        self.response = response
        self.period = period
        self.model = model
        self.offset = offset

    def evaluate(self, params, **user_args):
        from scipy.signal import fftconvolve
        model = self.model.evaluate(params, **user_args)
        n = len(model)
        shift = n % self.period
        n_periods = n // self.period + 1
        periodic_irf = np.roll(np.hstack(10*n_periods*[self.response[:self.period]]),
                               shift + self.offset)
        def pad(arr, n):
            m = arr.shape[0]
            if m < n:
                return np.hstack([arr, np.zeros(n-m)])
            else:
                return arr
        model = pad(model, n)
        a = fftconvolve(periodic_irf, model, 'same')
        return a[n_periods*self.period:n_periods*self.period+n]
    
    def parameters(self):
        return self.model.parameters()

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
polarization = True

irfs = [np.genfromtxt(irf, dtype=None, names='time,counts') for irf in args.irf]

if polarization:
    if len(irfs) != 2:
        raise RuntimeError('Expected two IRFs')
else:
    if len(irfs) != 1:
        raise RuntimeError('Expected one IRFs')

times = irfs[0]['time']
irfs = [irf['counts'] for irf in irfs]

# Determine the channel width (jiffy)
if args.jiffy is not None:
    jiffy = args.jiffy # in seconds
    jiffy_ps = jiffy / 1e-12
else:
    jiffy_ps = (times[1] - times[0]) # in picoseconds
    jiffy = jiffy_ps * 1e-12 # in seconds

# Determine the pulse repetition rate
if args.rep_rate is None:
    irf = irfs[0]
    middle = (np.max(irf) - np.min(irf)) / 2 + np.min(irf)
    idxs, = np.nonzero(np.logical_and(irf[:-1] < middle, middle < irf[1:]))
    a,b = sorted(idxs)[:2]
    per = b - a
    debug = False
    if debug:
        pl.plot(irf)
        pl.axhline(middle, c='k')
        pl.axvline(a, c='g')
        pl.axvline(b, c='g')
        pl.yscale('log')
        pl.show()
else:
    per = int(1 / args.rep_rate / jiffy) # period in ticks

print 'Period', per, 'bins'
print 'Channel width', jiffy_ps, 'ps'

n = args.periods * per
irfs = [irf[:n] for irf in irfs]

# Subtract background from and normalize IRF
def background_subtract(irf):
    bg = np.median(irf)
    print('IRF background = %1.2f' % bg)
    return irf - bg
irfs = map(background_subtract, irfs)

# Fix normalization of IRF
if polarization:
    # We need to take care to preserve the relative magnitude of the
    # perpendicular/parallel IRFs
    irfs[1] /= sum(irfs[1])
irfs[0] /= sum(irfs[0])

fit = Fit()

# Build decay model
rates = []
for i in range(args.components):
    tau = 1000 + 1000*i
    rate = fit.param('lambda%d' % i, initial=1/tau)
    rates.append(rate)

if polarization:
    r0 = fit.param('r0', initial=0.4)
    rate_rot = fit.param('lambda_rot', initial=1/1000)
    imbalance = fit.param('g', initial=1)

convolutions = []
for curve_idx,f in enumerate(args.corr):
    corr = np.genfromtxt(f)[:,1]
    assert len(corr) >= n
    corr = corr[:n] # FIXME?
    times = jiffy_ps * np.arange(corr.shape[0])
    weights = np.zeros_like(corr)
    weights[corr != 0] = 1 / np.sqrt(corr[corr != 0])
    name = f
    norm = np.sum(corr)

    if polarization:
        irf = irfs[curve_idx % 2]
        if curve_idx % 2 == 0:
            # Parallel channel
            rot_model = 1 + 2 * r0 * np.exp(-rate_rot * times)
            name += 'par'
            norm2 = 1
        else:
            # Perpendicular channel
            rot_model = 1 - r0 * np.exp(-rate_rot * times)
            name += 'perp'
            norm2 = imbalance
    else:
        rot_model = 1
        irf = irfs[0]

    decay_models = []
    for comp_idx, rate in enumerate(rates):
        amp = fit.param('c%d_amplitude%d' % (curve_idx, comp_idx),
                        initial=np.max(corr) / norm)
        decay_models.append(ExponentialModel(rate=rate, amplitude=amp))
    decay_model = sum(decay_models)

    convolved = ConvolvedModel(irf, per, decay_model * rot_model, offset=0)
    convolutions.append(convolved)
    model = norm2 * norm * convolved
    fit.add_curve(f, model, corr, weights=weights, t=times)

def fit_with_offset(offset):
    for m in convolutions:
        m.offset = offset
    return fit.fit()

offsets = [0] if args.no_offset else range(-5, 5)
fits = {i: fit_with_offset(i) for i in offsets}
print 'offsets', {i: fit.total_chi_sqr for i,fit in fits.items()}
offset,res = min(fits.items(), key=lambda (_, fit): fit.total_chi_sqr)
print 'optimal offset', offset

def print_params(p):
    if polarization:
        print '  g', p['g']
        print '  r0', p['r0']
        print '  tau_rot', 1/p['lambda_rot']

    for comp_idx in range(args.components):
        rate = p['lambda%d' % comp_idx]
        print '  Component %d' % comp_idx
        print '    tau', 1/rate

    for curve_idx,name in enumerate(args.corr):
        print '  Curve %s' % name
        for comp_idx in range(args.components):
            rate = p['lambda%d' % comp_idx]
            amp = p['c%d_amplitude%d' % (curve_idx, comp_idx)] / rate
            print '    amplitude%d' % comp_idx, amp

print
print 'Initial parameters'
print_params(res.initial.params)

print
print 'Fitted parameters'
print_params(res.params)


# Fix covariance
for comp_idx1 in range(args.components):
    for curve_idx1,_ in enumerate(args.corr):
        for comp_idx2 in range(args.components):
            for curve_idx2,_ in enumerate(args.corr):
                p1 = 'c%d_amplitude%d' % (curve_idx1, comp_idx1)
                p2 = 'c%d_amplitude%d' % (curve_idx2, comp_idx2)
                rate1 = res.params['lambda%d' % comp_idx1]
                rate2 = res.params['lambda%d' % comp_idx2]
                #res.covar[p1][p2] *= rate1 * rate2
    
print
print 'Reduced chi-squared'
for name, curve in res.curves.items():
    print '  %-15s     %1.3g' % (name, curve.reduced_chi_sqr)

print
print 'Standard error'
if res.stderr:
    for param, err in res.stderr.items():
        print '  %-15s     %1.2g' % (param, err)
else:
    print "  Failed to compute due to flat axis"

print
print 'Correlations (coefficients less than 0.2 omitted)'
if res.correl:
    correls = {(param1,param2): res.correl[param1][param2]
               for param1 in res.params.keys()
               for param2 in res.params.keys()
               if param1 < param2}
    for (p1,p2), c in sorted(correls.items(), key=lambda ((a,b),c): c, reverse=True):
        if abs(c) > 0.2:
            print '  %-15s / %-15s       %1.2f' % (p1, p2, c)
else:
    print "  Failed to compute due to flat axis"

plots = pl.subplot(211)
residuals = pl.subplot(212)
color_cycle = pl.rcParams['axes.color_cycle']
for curve_idx, name in enumerate(args.corr):
    cres = res.curves[name]
    #plots.plot(times, res.initial.curves[name].fit, '+', label='Initial')
    sym = '+' if curve_idx % 2 == 0 else 'x'
    color = color_cycle[curve_idx // 2 % len(color_cycle)]
    plots.plot(times, cres.curve.data, sym, label='Observed', color=color)
    plots.plot(times, cres.fit, label='Fit', color=color)
    residuals.plot(times, res.curves[name].residuals, color=color)

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
