#!/usr/bin/python

import numpy as np
from photon_tools.timetag_types import *
from alex import get_alex_bins
from matplotlib import pyplot as pl

clock = 1
binlength = 20
total_length = 1000
snippets = total_length / binlength


DD = 1
AD = 2
strobe_D = []
for i in np.arange(0, snippets):
    for j in np.arange(0, i%4+1):
        strobe_D.append(2*binlength*i + 2 + j)
    for j in np.arange(0, (5-i)%5+1):
        strobe_D.append(2*binlength*i + binlength + 2 + j)

strobe_D = np.core.records.fromarrays([np.uint64(strobe_D), np.uint8(strobe_D)], names='t,ch')

DA = 3
AA = 4
strobe_A = []
for i in np.arange(0, snippets):
    for j in np.arange(0, DA):
        strobe_A.append(2*binlength*i + 2 + j)
    for j in np.arange(0, AA):
        strobe_A.append(2*binlength*i + binlength + 2 + j)

strobe_A = np.core.records.fromarrays([np.uint64(strobe_A), np.uint8(strobe_A)], names='t,ch')

delta_D = np.empty(210, dtype=delta_event_dtype)
delta_D['start_t'] = binlength*np.arange(210)
delta_D['state'] = (np.arange(210)+1) % 2

delta_A = np.empty(210, dtype=delta_event_dtype)
delta_A['start_t'] = binlength*np.arange(210)
delta_A['state'] = (np.arange(210)) % 2

bins = get_alex_bins(strobe_D, strobe_A, delta_D, delta_A, binlength, clock)
# TODO: Check output

F_Dem_Dexc = bins.Dem_Dexc['count']
F_Dem_Aexc = bins.Dem_Aexc['count']
F_Aem_Dexc = bins.Aem_Dexc['count']
F_Aem_Aexc = bins.Aem_Aexc['count']

npts = 1000
pl.plot(F_Dem_Dexc[:npts], drawstyle='steps', label='Dem Dexc')
pl.plot(F_Dem_Aexc[:npts], drawstyle='steps', label='Dem Aexc')
pl.plot(F_Aem_Dexc[:npts], drawstyle='steps', label='Aem Dexc')
pl.plot(F_Aem_Aexc[:npts], drawstyle='steps', label='Aem Aexc')
pl.legend()
pl.show()
