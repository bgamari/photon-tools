#!/usr/bin/python

from photon_tools.timetag_types import *
from photon_tools.filter_photons import filter_by_spans
import numpy as np

strobe = np.empty(1e2, dtype=strobe_event_dtype)
strobe['t'] = np.arange(0,600,6)
strobe['chs'] = strobe['t']

delta_a = np.empty(15, dtype=delta_event_dtype)
delta_a['start_t'] = 10*np.arange(15)
delta_a['state'] = (np.arange(15)+1) % 2
#delta_a['start_t'][0] = 0
filtered_a = filter_by_spans(strobe, delta_a)
print(delta_a)
print(filtered_a)
print()

delta_d = np.empty(15, dtype=delta_event_dtype)
delta_d['start_t'] = 10*np.arange(15)
delta_d['state'] = np.arange(15) % 2
filtered_d = filter_by_spans(strobe, delta_d)

#print strobe
print(delta_d)
print(filtered_d)
