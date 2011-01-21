#!/usr/bin/python

import subprocess
import tempfile
import numpy as np
from numpy import any, all

def bin_data(data, bin_len):
        fi = tempfile.TemporaryFile('rw+b')
        di.tofile(fi)
        fi.seek(0)

        fo = tempfile.TemporaryFile('rw+b')
        subprocess.check_call(['./bin_photons', str(bin_len)], stdin=fi, stdout=fo)
        fo.seek(0)
        open('asdf', 'w').write(fo.read())
        fo.seek(0)
        dt = np.dtype([('time', 'u8'), ('counts','u2')])
        return np.fromfile('asdf', dtype=dt)

test_n = 1
def result(result, description):
        global test_n
        print 'ok' if result else 'not ok',
        print test_n, '-', description
        test_n += 1

print '1..3'

di = np.arange(1001, dtype='u8')
do = bin_data(di, 10)
result(all(do['counts'][0:99] == 10), 'basic binning')

di = np.array([0, 100, 110], dtype='u8')
do = bin_data(di, 10)
res = do['counts'][0] == 1 and all(do['counts'][1:9] == 0) and do['counts'][10] == 1
result(res, 'bin_photons zero bins test')

