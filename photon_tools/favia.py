import logging
import sys
from tempfile import NamedTemporaryFile
import subprocess
import numpy as np

#logging.basicConfig(level=logging.DEBUG)

dtype = np.dtype([('lag', 'f8'), ('loglag', 'f8'),
                  ('dot', 'f8'), ('dotnormed', 'f8'),
                  ('bar', 'f8')])

def read_favia(fname):
        return np.loadtxt(fname, dtype=dtype)
        
def acorr(x, jiffy=1./128e6, short_grain=1e-6, long_lag=1, verbose=False):
        return corr(x, x, jiffy, short_grain, long_lag, verbose)

def corr(x, y, jiffy=1./128e6, short_grain=1e-6, long_lag=1, verbose=False):
        """ Compute the correlation function of the datasets x and y. Jiffy,
        short_grain, and long_lag are given in seconds """
        fo = NamedTemporaryFile()
        fx = NamedTemporaryFile()
        fy = NamedTemporaryFile()
        x.tofile(fx.name)
        y.tofile(fy.name)

        args = ['favia',
                '--xfile=%s' % fx.name, '--yfile=%s' % fy.name,
                '--jiffy=%e' % jiffy,
                '--long_lag=%e' % long_lag,
                '--short_grain=%e' % short_grain]
        logging.debug(' '.join(args))
        stderr = sys.stderr if verbose else subprocess.PIPE
        p = subprocess.Popen(args, stdout=fo, stderr=stderr)
        if p.wait() != 0:
                if not verbose: print p.stderr.read()
                raise RuntimeError('Favia threw error: exit code %d' % p.returncode)

        return read_favia(fo.name)

