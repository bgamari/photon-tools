#!/usr/bin/python

import numpy as np
from photon_tools.timetag_types import *
from alex import get_alex_photons, get_alex_bins
from matplotlib import pyplot as pl

n_repeats = 40

dt_Dexc_Dem = 10*np.ones(100/10)
dt_Aexc_Aem = 20*np.ones(100/20)
dt_Aexc_Dem = 4*np.ones(100/4)
dt_Dexc_Aem = 5*np.ones(100/5)

strobe_Dem = np.cumsum(np.concatenate([dt_Dexc_Dem, dt_Aexc_Dem]*n_repeats))
strobe_Aem = np.cumsum(np.concatenate([dt_Dexc_Aem, dt_Aexc_Aem]*n_repeats))

strobe_D = np.core.records.fromarrays([np.uint64(strobe_Dem), np.uint8(strobe_Dem)], names='t,ch')
strobe_A = np.core.records.fromarrays([np.uint64(strobe_Aem), np.uint8(strobe_Aem)], names='t,ch')

delta_ts = np.arange(2*n_repeats)
delta_D = np.rec.array(zip(100*delta_ts, delta_ts % 2),
		       dtype=delta_event_dtype)
delta_A = np.rec.array(zip(100*delta_ts, (delta_ts+1) % 2),
		       dtype=delta_event_dtype)

def plot_test():
        pl.plot(strobe_A['t'], 6+np.ones_like(strobe_A['t']), 'ro', label='Acceptor strobe')
        pl.plot(strobe_D['t'], 5+np.ones_like(strobe_D['t']), 'b^', label='Donor strobe')
        pl.plot(delta_D['start_t'], 3+delta_D['state'], 'bo', label='Donor delta')
        pl.plot(delta_A['start_t'], 1+delta_A['state'], 'r^', label='Acceptor delta')
        for i in delta_ts:
                pl.axvline(x=100*i, alpha=0.4)
        pl.legend()
        pl.xlabel('Time')
        pl.xlim(0,1000)
        pl.ylim(0,9)
        pl.show()

photons = get_alex_photons(strobe_D, strobe_A, delta_D, delta_A, 0)
def iat(photons):
        """ Compute interarrival times """
        return photons['t'][1:] - photons['t'][:-1]
assert all(iat(photons.Dem_Dexc) == 4)
assert all(iat(photons.Aem_Aexc) == 5)
assert all(iat(photons.Aem_Dexc) == 20)
assert all(iat(photons.Dem_Aexc) == 10)

# The rest of this isn't being used at the moment
print 'good'
import sys
sys.exit(0)

bins = get_alex_bins(photons, 50)
F_Dem_Dexc = bins.Dem_Dexc['count']
F_Dem_Aexc = bins.Dem_Aexc['count']
F_Aem_Dexc = bins.Aem_Dexc['count']
F_Aem_Aexc = bins.Aem_Aexc['count']

