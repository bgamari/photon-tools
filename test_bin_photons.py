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
        dt = np.dtype([('time', 'u8'), ('counts','u2')])
        return np.fromfile(fo, dtype=dt)

test_n = 1
def result(result, description):
        global test_n
        print 'ok' if result else 'not ok',
        print test_n, '-', description
        test_n += 1

print '1..3'

di = np.arange(1001, dtype='u8')
do = bin_data(di, 10)
good = len(do) == 100 and all(do['counts'][0:100] == 10)
result(good, 'basic binning')

# Check that partial bins are thrown away
di = np.arange(1000, 1555, dtype='u8')
do = bin_data(di, 100)
good = len(do) == 5 \
        and all(do['counts'] == 100) \
        and all(do['time'] == np.arange(1000,1500,100))
result(good, 'partial bins')

di = np.array([0, 100, 110], dtype='u8')
do = bin_data(di, 10)
good = len(do) == 11 \
        and do['counts'][0] == 1 \
        and all(do['counts'][1:9] == 0) \
        and do['counts'][10] == 1
result(good, 'bin_photons zero bins test')

