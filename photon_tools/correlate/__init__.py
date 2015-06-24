from __future__ import division
from . import favia
from . import hphoton
import numpy as np

default_engine = 'favia'

def corr(x, y, jiffy, min_lag=1e-6, max_lag=1, fineness=8,
         engine=None, **kwargs):
    """
    Compute the correlation function of the photon timeseries ``x`` and
    ``y``. ``x`` and ``y`` should be arrays of unsigned integer
    timestamps with timestep `jiffy` (given in seconds). The
    correlation function will be computed for logarithmically-spaced
    lags starting from `short_grain` up to `long_lag` (given in
    seconds). The `fineness` parameter gives the number of lags
    computed per octave.

    :type x: array of integer timestamps
    :param x: Timeseries to convolve
    :type y: array of integer timestamps
    :param x: Timeseries to convolve (compute autocorrelation if ``None``)
    :type jiffy: ``float``, optional
    :param jiffy: The timestamp resolution
    :type min_grain: ``float``, optional
    :param min_grain: The shortest lag to compute
    :type max_lag: ``float``
    :param max_lag: The longest lag to compute
    :type fineness: ``int``, optional
    :param fineness: The number of lags to compute per octave
    :type engine: 'favia' or 'hphoton'
    :param engine: Which correlator to use to compute the correlation function.
    :param kwargs: Parameters passed to engine.
    """
    if engine is None:
        engine = default_engine

    if engine == 'favia':
        return favia.corr(x, y, jiffy, min_lag, max_lag, fineness, **kwargs)
    elif engine == 'hphoton':
        return hphoton.corr(x, y, jiffy, min_lag, max_lag, fineness, **kwargs)
    else:
        raise ValueError("Unknown correlator engine '%s'" % engine)

def autocorr(x, **kwargs):
    """
    Compute an auto-correlation function of the dataset ``x``. See
    :func:`corr` for a complete description of available options.
    """
    return corr(x, None, **kwargs)

def _split_at(timestamps, splits):
    """
    Split up an array into chunks.

        >>> _split_at(np.arange(5, 19), [11, 14])
        [[5,6,7,8,9,10], [11,12,13], [14,15,16,17,18]]

    :type timestamps: an :class:`array`s of timestamps
    :type splits: iterable
    """
    # Since numpy doesn't yet offer a `find` function, we need
    # take effort to avoid doing repeated work
    chunks = []
    xs = timestamps
    for i, upper in enumerate(splits):
        take, = np.nonzero(xs >= upper)
        if len(take) == 0:
            chunks.append(xs)
            xs = xs[:0]
        else:
            idx = take[0]
            chunks.append(xs[:idx])
            xs = xs[idx:]

    chunks.append(xs)
    return chunks

def corr_chunks(x, y, n=10, cross_chunks=False, anomaly_thresh=None, **kwargs):
    """
    Compute the cross-correlation between two photon timeseries ``x``
    and ``y``, computing the variance by splitting the series into
    ``n`` chunks.

    :type n: ``int``
    :param n: The number of chunks to use for computation of the variance.
    :param kwargs: Keyword arguments to be passed to :func:`corr`
    :type cross_chunks: :class:`bool`
    :param cross_chunks: Use cross-correlations between non-cooccurrant chunks.
    :type anomaly_thresh: (optional) :class:`float`
    :param anomaly_thresh: Anomaly normalized-log-likelihood threshold. See :func:`anomaly_thresh`
    :returns: tuple of ``(Gmean, corrs)`` where ``Gmean`` is a record array with
      fields ``lag``, ``G``, and ``var`` and ``corrs`` is a :class:`array` of
      shape ``(n,nlags)`` containing the correlation functions of the individual
      chunks.
    """
    max_t = max(x[-1], y[-1])
    min_t = min(x[0], y[-1])
    dt = (max_t - min_t) / (n+1)
    splits = np.arange(1, n) * dt + min_t

    x_chunks = _split_at(x, splits)
    y_chunks = _split_at(y, splits)

    if cross_chunks:
        pairs = [(x, y) for x in x_chunks for y in y_chunks]
    else:
        pairs = zip(x_chunks, y_chunks)

    # g.shape == (Nlags,)
    g = corr(x, y, **kwargs)
    lags = g['lag']
    g = g['G']
    # corrs.shape == (len(pairs), Nlags)
    corrs = np.vstack( corr(xc, yc, **kwargs)['G'] for (xc,yc) in pairs )
    corrs = (corrs-1) * ((g-1).sum() / (corrs-1).sum(axis=1))[:,np.newaxis] + 1

    if anomaly_thresh is not None:
        likelihoods = anomaly_likelihood(corrs) / corrs.shape[1]
        take = likelihoods > anomaly_thresh
        if np.count_nonzero(take) == 0:
            raise ValueError('No chunks deemed non-anomalous')
        corrs = corrs[take, :]

        g = np.mean(corrs, axis=0)

    var = np.var(corrs, axis=0) / n

    return (np.rec.fromarrays([lags, g, var], names='lag,G,var'), corrs)

def anomaly_likelihood(xs):
    """
    Evaluate the likelihood of each sample under a Gaussian model induced by the
    others.

    :type xs: :class:`array` of shape ``(Nsamples, Nfeatures)``
    :param xs: samples
    :rtype: :class:`array` of shape ``(Nsamples,)``
    :returns: The log likelihood
    """
    likelihoods = []
    for i in range(xs.shape[0]):
        take = np.ones(xs.shape[0], dtype=bool)
        take[i] = False
        subsample = xs[take,:]
        mean = np.mean(subsample, axis=0)
        var = np.var(subsample, axis=0)

        test = xs[i,:]
        likelihood = -(test - mean)**2 / 2 / var**2 - 0.5 * np.log(2*np.pi*var)
        likelihoods.append(np.sum(likelihood))

    return np.array(likelihoods)
