from __future__ import division
import numpy as np
from matplotlib import pyplot as pl
import squmfit
from squmfit import Fit, model, Argument
from collections import namedtuple

def make_map(type):
    def map(self, f):
        return type(*[f(x) for x in self])
    setattr(type, "map", map)

class Aniso(object):
    """ A pair of anisotropy detection channels """
    def __init__(self, par, perp):
        self.par = par
        self.perp = perp

    def map(self, f):
        return Aniso(f(self.par), f(self.perp))

class FitSet(object):
    """ A pair of anisotropy decay curves along with their IRF """
    def __init__(self, name, irf, decay):
        """
        :type name: str
        :type irf: :class:`Aniso`
        :param irf: normalized IRF histogram
        :type decay: :class:`Aniso`
        """
        self.name = name
        self.irf = irf
        self.decay = decay

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
    #from scipy.signal import fftconvolve
    #a = fftconvolve(model, response, 'valid')
    from scipy.fftpack import fft, ifft
    a = ifft(fft(response) * fft(model)).real
    debug = False
    if debug:
        pl.plot(response)
        pl.plot(model)
        pl.plot(a)
        #pl.yscale('log')
        pl.show()
    return a

def estimate_rep_rate(irf, debug=False):
    """
    Estimate excitation repetition rate (measured in bins) from IRF histogram.

    :param irf: IRF histogram counts
    :type irf: array of shape ``(N,)``
    :rtype: int
    """
    middle = (np.max(irf) - np.min(irf)) / 2 + np.min(irf)
    idxs, = np.nonzero(np.logical_and(irf[:-1] < middle, middle < irf[1:]))
    a,b = sorted(idxs)[:2]
    period = b - a
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

    :type irfs: :class:`Aniso` of histograms
    :rtype: :class:`Aniso` of histograms
    """
    def background_subtract(irf):
        bg = np.median(irf)
        print('IRF background = %1.2f' % bg)
        return irf - bg

    irfs = irfs.map(background_subtract)

    # Fix normalization of IRF
    return irfs.map(lambda x: x / sum(x))

def fit(corrs, jiffy_ps, exc_period, n_components, periods=1, **kwargs):
    """
    :type corrs: list of :class:`FitSet`s
    :param corrs: Fluorescence histograms along with their IRFs
    :type jiffy_ps: int
    :param jiffy_ps: Bin width in picoseconds
    :type exc_period: int
    :param exc_period: Excitation period measured in bins
    :type n_components: int
    :param n_components: Number of exponential decay components to fit against
    :type periods: int
    :param periods: Number of periods to fit against
    :param kwargs: Other keyword arguments passed to `analyze`
    :rtype: tuple of (:class:`FitResults`, :class:`FitResults`, :class:`CurveDesc`)
    :returns:
      1. :class:`FitResult` of fit with period held fixed at initial value
      2. :class:`FitResult` of fit with period free
      3. :class:`CurveDesc` of fit with period free
    """
    # Run the fit first to get the parameters roughly correct, then
    # then infer the period
    res1,desc1 = analyze(corrs, exc_period, n_components, jiffy_ps, **kwargs)
    res2,desc2 = analyze(corrs, exc_period, n_components, jiffy_ps,
                   free_period=True, params0=res1.params, **kwargs)
    return res1, res2, desc2

# curve parameters
CurveDesc = namedtuple('CurveDesc', 'amps exc_leakage tau_rot')
# model parameters
ModelDesc = namedtuple('ModelDesc', 'fluor_rates period offset_par offset_perp r0 imbalance tau_rot curves')

def analyze(corrs, exc_period, n_components, jiffy_ps,
            params0=None, free_period=False,
            exc_leakage=False, imbalance=None, indep_aniso=False,
            no_offset=False, fix_lifetimes=[]):
    """
    Fit a set of anisotropy data with the given IRF and model

    :type corrs: :class:`FitSet`
    :type exc_period: `float`
    :param exc_period: the approximate excitation period in bins
    :type n_components: `int`
    :param n_components: the number of exponential components to fit against
    :type jiffy_ps: `int`
    :param jiffy_ps: the channel width in picoseconds
    :param params0: The initial parameters to fit with
    :type free_period: `bool`
    :param free_period: Whether to fit the excitation period
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
    :type fix_lifetimes: list of floats, times in picoseconds
    :param fix_lifetimes: Fix the lifetimes of some fluorescence decay components
    """
    fit = Fit()
    lag = Argument('t')

    offset_par = 0 if no_offset else fit.param('offset-par', 0)
    offset_perp = 0 if no_offset else fit.param('offset-perp', 0)
    period = fit.param('period', exc_period)

    # Build decay model
    assert len(fix_lifetimes) <= n_components
    rates = []
    for i in range(n_components):
        if i < len(fix_lifetimes):
            rate = 1 / fix_lifetimes[i]
        else:
            tau = 1000 + 1000*i
            rate = fit.param('lambda%d' % i, initial=1/tau)
        rates.append(rate)

    # Parameters for anisotropy model
    r0 = fit.param('r0', initial=0.4)
    if imbalance is None:
        imbalance = fit.param('g', initial=1)
    if not indep_aniso:
        tau_rot = fit.param('tau_rot', initial=500)

    curve_descs = []
    for pair_idx,corr in enumerate(corrs):
        if indep_aniso:
            tau_rot = fit.param('%s_tau_rot' % pair.name, initial=500)

        # generate fluorescence decay model
        n = len(corr.irf.par) # FIXME?
        assert len(corr.decay.par) >= n
        assert len(corr.decay.perp) >= n
        par = corr.decay.par
        perp = corr.decay.perp

        decay_models = []
        initial_amp = np.max(par) / np.sum(par)
        amplitudes = []
        for comp_idx, rate in enumerate(rates):
            amp = fit.param('%s_amplitude%d' % (corr.name, comp_idx), initial=initial_amp)
            amplitudes.append(amp / rates[comp_idx])
            decay_models.append(exponential(t=lag, rate=rate, amplitude=amp))
        decay_model = sum(decay_models)

        # IRF leakage
        leakage = 0
        if exc_leakage:
            leakage = fit.param('%s_leakage' % corr.name, initial=0.1)

        def add_curve(corr, name, rot_model, offset, norm, irf):
            times = jiffy_ps * np.arange(10*corr.shape[0])

            # generate weights
            weights = np.zeros_like(corr, dtype='f')
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
                  name = '%s_par' % corr.name,
                  irf = corr.irf.par)
        add_curve(corr = perp,
                  rot_model = lambda times: 1 - r0 * np.exp(-1 / tau_rot * times),
                  norm = imbalance,
                  offset = offset_perp,
                  name = '%s_perp' % corr.name,
                  irf = corr.irf.perp)

        curve_descs.append(CurveDesc(
            amps = amplitudes,
            tau_rot = tau_rot,
            exc_leakage = leakage,
        ))

    model_desc = ModelDesc(
        fluor_rates = rates,
        period = period,
        offset_par=offset_par,
        offset_perp = offset_perp,
        r0 = r0,
        imbalance = imbalance,
        tau_rot = tau_rot,
        curves = curve_descs,
    )
    res = fit.fit(params0)
    return (res, model_desc)

def plot(fig, corrs, jiffy_ps, result, sep_resid=False, opacity=0.4):
    """
    Plot a fit produced by :func:`analyze`

    :type fig: :class:`pl.Figure`
    :type corrs: :class:`Aniso` of :class:`file`s
    :type jiffy_ps: float
    :param jiffy_ps: Jiffy in picoseconds
    :param result: A :class:`FitResult` from :func:`analyze`
    :type sep_resid: `bool`
    :param sep_resid: Plot residuals on per-pair axes.
    """
    import matplotlib.gridspec as gridspec
    if sep_resid:
        ncurves = len(corrs)
        gs = gridspec.GridSpec(1 + ncurves, 2,
                               width_ratios=[3,1],
                               height_ratios=[4] + [1]*ncurves)
    else:
        gs = gridspec.GridSpec(2, 2, width_ratios=[3,1], height_ratios=[3,1])

    plots = pl.subplot(gs[0, 0])
    legend = gs[0:1, 1]
    if sep_resid:
        residuals = {pair_idx: pl.subplot(gs[1+pair_idx, 0]) for pair_idx,_ in enumerate(corrs)}
    else:
        resid = pl.subplot(gs[1, 0])
        residuals = {pair_idx: resid for pair_idx,_ in enumerate(corrs)}

    color_cycle = pl.rcParams['axes.color_cycle']
    for pair_idx, aniso in enumerate(corrs):
        for ch in ['par', 'perp']:
            cres = result.curves['%s_%s' % (aniso.name, ch)]
            color = color_cycle[pair_idx % len(color_cycle)]
            times = jiffy_ps / 1000 * np.arange(cres.fit.shape[0])
            sym = '+' if ch == 'par' else 'x'
            label = aniso.name if ch == 'par' else None
            kwargs = {'color': color, 'markersize': 1.5}
            plots.plot(times, cres.curve.data, sym, alpha=opacity/2, **kwargs)
            plots.plot(times, cres.fit, label=label, alpha=opacity, **kwargs)
            residuals[pair_idx].plot(times, cres.residuals, sym, alpha=opacity/3, **kwargs)

    plots.set_ylabel('number of occurences')
    for k,axes in residuals.items():
        axes.axhline(0, color='k')
        pl.setp(axes.get_xticklabels(), visible=False)
        if sep_resid:
            axes.locator_params(axis='y', nbins=2)

    residuals.values()[len(corrs)//2].set_ylabel('residual')
    residuals.values()[len(corrs)-1].set_xlabel('time (ns)')
    pl.setp(residuals.values()[len(corrs)-1].get_xticklabels(), visible=True)
    plots.set_yscale('log')
    pl.setp(plots.get_xticklabels(), visible=False)

    import matplotlib.lines as mlines
    handles, labels = plots.get_legend_handles_labels()
    handles = [
        mlines.Line2D([], [], color='k', label='fit'),
        mlines.Line2D([], [], color='k', marker='+', label='parallel'),
        mlines.Line2D([], [], color='k', marker='x', label='perpendicular'),
        mlines.Line2D([], [], alpha=0, label=''),  # Spacer
    ] + handles
    bbox = legend.get_position(fig)
    plots.legend(handles=handles, loc='upper left',
                 bbox_to_anchor=bbox,
                 bbox_transform=fig.transFigure,
                 mode='expand', fontsize='small', ncol=1, frameon=False)
