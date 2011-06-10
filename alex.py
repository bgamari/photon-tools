#!/usr/bin/python

import argparse
import numpy as np
from photon_tools.filter_photons import filter_by_spans
from photon_tools.timetag_parse import get_strobe_events, get_delta_events
from photon_tools.bin_photons import bin_photons

parser = argparse.ArgumentParser()
parser.add_argument('file', type=argparse.FileType('r'),
                    help='Timetag file')
parser.add_argument('-s', '--start-offset', type=float,
                    help='Time offset of valid data after delta channel goes high (seconds)', default=1e-6)
parser.add_argument('-w', '--bin-width', type=float,
                    help='Bin width (seconds)', default=1e-3)
args = parser.parse_args()

start_exc_offset = args.start_offset
bin_width = 1e-3
strobe_clock = delta_clock = 128e6
f = args.file.name

Dem = get_strobe_events(f, 0x1)
Aem = get_strobe_events(f, 0x2)
delta_D = get_delta_events(f, 0)
delta_A = get_delta_events(f, 1)

def shifted_deltas(deltas, state, off):
        """ Add an offset to the start times of delta records with the given state """
        ret = np.array(deltas)
        taken = ret['state'] == state
        ret['start_t'][taken] += off*delta_clock
        return ret

Dem_Dexc = filter_by_spans(strobe_D, shifted_deltas(strobe_D, True, start_exc_offset))
Dem_Aexc = filter_by_spans(strobe_D, shifted_deltas(strobe_A, True, start_exc_offset))
Aem_Dexc = filter_by_spans(strobe_A, shifted_deltas(strobe_D, True, start_exc_offset))
Aem_Aexc = filter_by_spans(strobe_A, shifted_deltas(strobe_A, True, start_exc_offset))

Dem_Dexc_bins = bin_photons(Dem_Dexc['t'], bin_width*strobe_clock)
Dem_Aexc_bins = bin_photons(Dem_Aexc['t'], bin_width*strobe_clock)
Aem_Dexc_bins = bin_photons(Aem_Dexc['t'], bin_width*strobe_clock)
Aem_Aexc_bins = bin_photons(Aem_Aexc['t'], bin_width*strobe_clock)

F_D_Dem_Dexc = Dem_Dexc_bins['count']
F_fret = Aem_Dexc_bins['count'] - Lk - Dir
F_A_Aem_Aexc = Aem_Aexc_bins['count']

