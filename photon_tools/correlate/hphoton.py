import logging
import sys
from tempfile import NamedTemporaryFile
import subprocess
import numpy as np

"""
Python interface to hphoton-correlate utility for computing
fluorescence correlation functions.

See <http://github.com/bgamari/hphoton>
"""

#logging.basicConfig(level=logging.DEBUG)
keep = False # For debugging

dtype = np.dtype([('lag', 'f8'), ('G', 'f8'), ('var', 'f8')])

class CorrelateError(RuntimeError):
    def __init__(self, exit_code, error):
        RuntimeError.__init__(self, 'hphoton correlate threw an error: exit code %d\n\n%s' % (exit_code, error))
        self.exit_code = exit_code
        self.error = error

def read_correlate(fname):
    return np.loadtxt(fname, dtype=dtype)

def corr(x, y, jiffy=1./128e6, min_lag=1e-6, max_lag=1, fineness=8, verbose=False):
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
    :param y: Timeseries to convolve (compute autocorrelation is ``None``)
    :type jiffy: ``float``, optional
    :param jiffy: The timestamp resolution
    :type min_lag: ``float``, optional
    :param min_lag: The shortest lag to compute
    :type max_lag: ``float``
    :param max_lag: The longest lag to compute
    :type fineness: ``int``, optional
    :param fineness: The number of lags to compute per octave
    :type verbose: ``bool``, optional
    :param verbose: Enable verbose output from ``favia`` correlator
    """
    fo = NamedTemporaryFile(delete=not keep)
    fx = NamedTemporaryFile(delete=not keep, suffix='.raw')
    x.tofile(fx.name)
    if y is not None:
        fy = NamedTemporaryFile(delete=not keep, suffix='.raw')
        y.tofile(fy.name)

    args = [
        'correlate',
        '-x'+fx.name,
        '-y'+fy.name if y is not None else '',
        '--jiffy=%e' % jiffy,
        '--max-lag=%e' % max_lag,
        '--min-lag=%e' % min_lag,
        '--nbins=%d' % fineness
    ]
    logging.debug(' '.join(args))
    stderr = sys.stderr if verbose else subprocess.PIPE
    p = subprocess.Popen(args, stdout=fo, stderr=stderr)
    error = p.stderr.read() if not verbose else None
    if p.wait() != 0:
        raise CorrelateError(p.returncode, error)

    return read_correlate(fo.name)
