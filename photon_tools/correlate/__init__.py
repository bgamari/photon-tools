import favia
import hphoton
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
    :param x: Timeseries to convolve
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
    :function:`corr` for a complete description of available options.
    """
    return corr(x, x, **kwargs)

def _split_chunks(x, n):
    l = len(x) / n
    return [ x[i*l:(i+1)*l] - x[i*l] for i in range(n) ]

def corr_chunks(x, y, n=10, **kwargs):
    """
    Compute the cross-correlation between two photon timeseries ``x``
    and ``y``, computing the variance by splitting the series into
    ``n`` chunks.

    :type n: ``int``
    :param n: The number of chunks to use for computation of the variance.
    :param kwargs: Keyword arguments to be passed to :func:`corr`
    """
    x_chunks = _split_chunks(x, n)
    y_chunks = _split_chunks(y, n)
    corrs = np.vstack( corr(xc, yc, **kwargs) for (xc,yc) in zip(x_chunks,y_chunks) )
    g = corr(x, y, **kwargs)['G']
    var = np.var(corrs['G'], axis=0) / n
    return (np.rec.fromarrays([corrs[0]['lag'], g, var], names='lag,G,var'), corrs['G'])
