import logging
import sys
from tempfile import NamedTemporaryFile
import subprocess
import numpy as np

"""
Python interface to favia utility for computing fluorescence
correlation functions.
"""

#logging.basicConfig(level=logging.DEBUG)
keep = False # For debugging

raw_dtype = np.dtype([('lag', 'f8'), ('loglag', 'f8'),
                      ('dot', 'f8'), ('dotnormed', 'f8'),
                      ('var', 'f8')])

class FaviaError(RuntimeError):
    def __init__(self, exit_code, error):
        RuntimeError.__init__(self, 'Favia threw an error: exit code %d\n\n%s' % (exit_code, error))
        self.exit_code = exit_code
        self.error = error

def read_favia_raw(fname):
    return np.loadtxt(fname, dtype=raw_dtype)

def read_favia(fname):
    c = read_favia_raw(fname)
    return np.rec.fromarrays([c['lag'], c['dotnormed'], c['var']],
                             names='lag,G,var')

def acorr(x, **kwargs):
    """
    Compute an auto-correlation function of the datasets x. See
    :function:`corr` for a complete description of available options.
    """
    return corr(x, x, **kwargs)

def corr(x, y, jiffy=1./128e6, short_grain=1e-6, long_lag=1, fineness=8, verbose=False):
    """
    Compute the correlation function of the datasets `x` and
    `y`. `x` and `y` should be unsigned integer timestamps with
    timestep `jiffy` (given in seconds). The correlation function
    will be computed for logarithmically-spaced lags starting from
    `short_grain` up to `long_lag` (given in seconds). The
    `fineness` parameter gives the number of lags computed per
    octave.

    :type x: array of integer timestamps
    :param x: Timeseries to convolve
    :type y: array of integer timestamps
    :param x: Timeseries to convolve
    :type jiffy: ``float``, optional
    :param jiffy: The timestamp resolution
    :type short_grain: ``float``, optional
    :param short_grain: The shortest lag to compute
    :type long_lag: ``float``
    :param long_lag: The longest lag to compute
    :type fineness: ``int``, optional
    :param fineness: The number of lags to compute per octave
    :type verbose: ``bool``, optional
    :param verbose: Enable verbose output from ``favia`` correlator
    """
    fo = NamedTemporaryFile(delete=not keep)
    fx = NamedTemporaryFile(delete=not keep)
    fy = NamedTemporaryFile(delete=not keep)
    x.tofile(fx.name)
    y.tofile(fy.name)

    args = [
            'favia',
            '--xfile=%s' % fx.name, '--yfile=%s' % fy.name,
            '--jiffy=%e' % jiffy,
            '--long_lag=%e' % long_lag,
            '--short_grain=%e' % short_grain,
            '--fineness=%d' % fineness
    ]
    logging.debug(' '.join(args))
    stderr = sys.stderr if verbose else subprocess.PIPE
    p = subprocess.Popen(args, stdout=fo, stderr=stderr)
    error = p.stderr.read() if not verbose else None
    if p.wait() != 0:
        raise FaviaError(p.returncode, error)

    return read_favia(fo.name)
