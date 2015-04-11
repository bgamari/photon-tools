from __future__ import division
import numpy as np
import squmfit
from squmfit import Fit, model, Argument
from collections import namedtuple

class Aniso(object):
    """ Anisotropy detection channels """
    def __init__(self, par, perp):
        self.par = par
        self.perp = perp

    def map(self, f):
        return Aniso(f(self.par), f(self.perp))

@model
def exponential(t, rate, amplitude):
    """ Note that this should contain a prefactor of `rate` that has
    been omitted to minimize parameter covariance """
    return amplitude * np.exp(-t * rate)

@model
def interpolate_irf(response, period, offset, periods=1):
    """
    An instrument response function which is interpolated and
    (possibly) shifted in time.

    Evaluates to an array of the same shape as `response` times `periods`.

    :type response: array of shape (N,)
    :param response: Instrument response
    :param period: The period of the instrument response
    :param offset: The temporal offset between the response and
        the measured curve (in bins)
    :type periods: int
    :param periods: The number of periods to produce
    """

    from scipy.interpolate import interp1d

    shift = len(response) % period
    periodic_irf_f = interp1d(np.arange(len(response)), response,
                              assume_sorted=True)
    ts = np.arange(3 * periods * len(response)) + shift + offset
    ts %= period
    # sometimes numerical error gives us values slightly outside
    # the interpolation range, be sure to clamp these
    ts = np.clip(ts, a_min=0, a_max = len(response) - 1)
    interpolated_irf = periodic_irf_f(ts)
    return interpolated_irf[:periods * len(response)]

@model
def convolved_model(response, model):
    """
    Convolve a model with a response.

    :type response: :class:`InterpolatedIrf`
    :param response: The instrument response
    :type model: :class:`squmfit.Expr`
    :param model: The model we are trying to evaluate
    """
    from scipy.signal import fftconvolve
    a = fftconvolve(response, model, 'same')
    return a

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
    return irfs.map(lambda x: x / sum(x))

def fit(irfs, corrs, jiffy_ps, exc_period, n_components, periods=1, **kwargs):
    """
    :type irfs: :class:`Aniso` of arrays
    :param irfs: Normalized IRF histograms
    :type corrs: list of :class:`Aniso`s
    :param corrs: Fluorescence histograms
    :type jiffy_ps: int
    :param jiffy_ps: Bin width in picoseconds
    :type exc_period: int
    :param exc_period: Excitation period measured in bins
    :type n_components: int
    :param n_components: Number of exponential decay components to fit against
    :type periods: int
    :param periods: Number of periods to fit against
    :param kwargs: Other keyword arguments passed to `analyze`
    :rtype: tuple of two :class:`FitResults`
    """
    jiffy = jiffy_ps * 1e-12
    n = periods * exc_period
    irfs = irfs.map(lambda x: x[:n])

    # Run the fit first to get the parameters roughly correct, then
    # then infer the period
    res1 = analyze(irfs, corrs, exc_period, n_components, jiffy_ps, **kwargs)
    res2 = analyze(irfs, corrs, exc_period, n_components, jiffy_ps,
                   free_period=True, params0=res1.params, **kwargs)
    return res1, res2

def analyze(irfs, corrs, exc_period, n_components, jiffy_ps,
            params0=None, free_period=False, trim_end=0,
            exc_leakage=False, imbalance=None, indep_aniso=False,
            no_offset=False):
    """
    Fit a set of anisotropy data with the given IRF and model

    :type irfs: :class:`Aniso` of arrays
    :param irfs: Normalized IRF histograms
    :type corrs: :class:`Aniso` of :class:`file`s
    :type exc_period: `float`
    :param exc_period: the approximate excitation period in bins
    :type n_components: `int`
    :param n_components: the number of exponential components to fit against
    :type jiffy_ps: `int`
    :param jiffy_ps: the channel width in picoseconds
    :param params0: The initial parameters to fit with
    :type free_period: `bool`
    :param free_period: Whether to fit the excitation period
    :type trim_end: `int`
    :param trim_end: Number of bins to drop off of the end of the correlations.
       This is to work around the spurious bins with low counts near the end of
       histograms produced by the Picoharp 300.
    :type exc_leakage: `bool`
    :param exc_leakage: Whether to include leakage of the excitation (modelled by
       the IRF) in the detection model.
    :type imbalance: `float` or None
    :param imbalance: Fix detector imbalance factor g
    :type indep_aniso: `bool`
    :param indep_aniso: Whether to fit a different rotational
       coherence time for each curve.
    :type no_offset: `bool`
    :param no_offset: Disable fitting of temporal offset between IRF and fluorescence curve.
    """
    fit = Fit()
    lag = Argument('t')

    offset_par = 0 if no_offset else fit.param('offset-par', 0)
    offset_perp = 0 if no_offset else fit.param('offset-perp', 0)
    period = exc_period if False else fit.param('period', exc_period)

    # Build decay model
    rates = []
    for i in range(n_components):
        tau = 1000 + 1000*i
        rate = fit.param('lambda%d' % i, initial=1/tau)
        rates.append(rate)

    # Parameters for anisotropy model
    r0 = fit.param('r0', initial=0.4)
    if imbalance is None:
        imbalance = fit.param('g', initial=1)
    if not indep_aniso:
        tau_rot = fit.param('tau_rot', initial=1000)

    for pair_idx,corr in enumerate(corrs):
        if indep_aniso:
            tau_rot = fit.param('tau_rot%d' % pair_idx, initial=1000)

        def read_histogram(path):
            # read histogram
            corr = np.genfromtxt(path)[:,1]
            n = len(irfs.par) # FIXME?
            assert len(corr) >= n
            corr = corr[:n - trim_end] # FIXME?
            return corr

        # generate fluorescence decay model
        par = read_histogram(corr.par.name)
        perp = read_histogram(corr.perp.name)
        assert len(par) == len(perp)

        decay_models = []
        initial_amp = np.max(par) / np.sum(par)
        for comp_idx, rate in enumerate(rates):
            amp = fit.param('c%d_amplitude%d' % (pair_idx, comp_idx), initial=initial_amp)
            decay_models.append(exponential(t=lag, rate=rate, amplitude=amp))
        decay_model = sum(decay_models)

        # IRF leakage
        leakage = 0
        if exc_leakage:
            leakage = fit.param('c%d_leakage' % pair_idx, initial=0.1)

        def add_curve(corr, name, rot_model, offset, norm, irf):
            times = jiffy_ps * np.arange(10*corr.shape[0])

            # generate weights
            weights = np.zeros_like(corr)
            weights[corr != 0] = 1 / np.sqrt(corr[corr != 0])

            # generate model
            iirf = interpolate_irf(response=irf,
                                   period=period if free_period else exc_period,
                                   offset=offset,
                                   periods=10)
            convolved = convolved_model(response=iirf, model=decay_model * rot_model(times))
            model = norm * convolved
            if leakage != 0:
                model = model + leakage * iirf.map(lambda a: a[:len(corr)])

            fit.add_curve(name, np.sum(corr) * model[:len(corr)], corr, weights=weights, t=times)

        add_curve(corr = par,
                  rot_model = lambda times: 1 + 2 * r0 * np.exp(-1 / tau_rot * times),
                  norm = 1,
                  offset = offset_par,
                  name = corr.par.name+'_par',
                  irf = irfs.par)
        add_curve(corr = perp,
                  rot_model = lambda times: 1 - r0 * np.exp(-1 / tau_rot * times),
                  norm = imbalance,
                  offset = offset_perp,
                  name = corr.perp.name+'_perp',
                  irf = irfs.perp)

    return fit.fit(params0)
