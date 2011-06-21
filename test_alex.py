#!/usr/bin/python

import numpy as np
from photon_tools.timetag_types import *
from alex import get_alex_photons, get_alex_bins
from matplotlib import pyplot as pl

n_repeats = 4

dt_Dexc_Dem = 10*np.ones(100/10)
dt_Aexc_Aem = 20*np.ones(100/20)
dt_Aexc_Dem = 2*np.ones(100/2)
dt_Dexc_Aem = 4*np.ones(100/4)

strobe_Dem = np.cumsum(np.concatenate([dt_Dexc_Dem, dt_Aexc_Dem]*n_repeats))
strobe_Aem = np.cumsum(np.concatenate([dt_Dexc_Aem, dt_Aexc_Aem]*n_repeats))

strobe_D = np.core.records.fromarrays([np.uint64(strobe_Dem), np.uint8(strobe_Dem)], names='t,ch')
strobe_A = np.core.records.fromarrays([np.uint64(strobe_Aem), np.uint8(strobe_Aem)], names='t,ch')

delta_ts = np.arange(2*n_repeats)
delta_D = np.rec.array(zip(100*delta_ts, delta_ts % 2),
		       dtype=delta_event_dtype)
delta_A = np.rec.array(zip(100*delta_ts, (delta_ts+1) % 2),
		       dtype=delta_event_dtype)

photons = get_alex_photons(strobe_D, strobe_A, delta_D, delta_A, 0)
bins = get_alex_bins(photons, 50)
# TODO: Check output

F_Dem_Dexc = bins.Dem_Dexc['count']
F_Dem_Aexc = bins.Dem_Aexc['count']
F_Aem_Dexc = bins.Aem_Dexc['count']
F_Aem_Aexc = bins.Aem_Aexc['count']

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

def plot_iat(photons, label):
	t = photons['t']
	print label, t, t[1:] - t[:-1], t.dtype
	#pl.plot(t[1:] - t[:-1], label=label)
	pl.plot(photons['t'], '+', label=label)
plot_iat(photons.Dem_Dexc, 'Dem_Dexc')
plot_iat(photons.Dem_Aexc, 'Dem_Aexc')
plot_iat(photons.Aem_Dexc, 'Aem_Dexc')
plot_iat(photons.Aem_Aexc, 'Aem_Aexc')
#pl.ylim(0, 10)
pl.ylabel('Photon IAT')
pl.xlabel('Photon Number')
pl.legend()
pl.show()

print F_Dem_Dexc
pl.plot(F_Dem_Dexc, drawstyle='steps', label='Dem Dexc')
pl.plot(F_Dem_Aexc, drawstyle='steps', label='Dem Aexc')
pl.plot(F_Aem_Dexc, drawstyle='steps', label='Aem Dexc')
pl.plot(F_Aem_Aexc, drawstyle='steps', label='Aem Aexc')
pl.legend()
pl.show()
