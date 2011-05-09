from tempfile import NamedTemporaryFile
import subprocess
import numpy as np

dtype = np.dtype([('lag', 'f8'), ('loglag', 'f8'),
                  ('dot', 'f8'), ('dotnormed', 'f8'),
                  ('bar', 'f8')])

def acorr(x, clockrate=160e6, short_grain=1e-6, long_lag=1):
        return corr(x, x, clockrate, short_grain, long_lag)

def corr(x, y, clockrate=160e6, short_grain=1e-6, long_lag=1):
        fx = NamedTemporaryFile()
        fy = NamedTemporaryFile()
        np.tofile(fx, x)
        np.tofile(fy, y)
        args = ['favia',
                '--xfile=%s ' % fx.name, '--yfile=%s' % fy.name,
                '--jiffy=%f' % 1./clockrate,
                '--long_lag=%u' % long_lag,
                '--short_grain=%u' % short_grain]
        fo = NamedTemporaryFile()
        p = subprocess.Popen(args, stdout=fo)
        out = p.communicate()
        if p.returncode != 0:
                raise RuntimeError('Favia threw error')

        return np.loadtxt(fo, dtype=dtype)

