#!/usr/bin/python

import argparse
import numpy as np
from photon_tools.filter_photons import filter_by_spans
from photon_tools.timetag_parse import get_strobe_events, get_delta_events
from photon_tools.bin_photons import bin_photons
from photon_tools.timetag_types import *
from matplotlib import pyplot as pl
from numpy import mean, std, amin, amax, logical_and

parser = argparse.ArgumentParser()
parser.add_argument('file', type=argparse.FileType('r'),
                    help='Timetag file')
parser.add_argument('-s', '--start-offset', type=float,
                    help='Time offset of valid data after delta channel goes high (seconds)', default=1e-6)
parser.add_argument('-w', '--bin-width', type=float,
                    help='Bin width (seconds)', default=1e-3)
parser.add_argument('-o', '--output', metavar='FILE', help='Output File Name')

args = parser.parse_args()

start_exc_offset = args.start_offset
bin_width = args.bin_width
strobe_clock = delta_clock = 128e6
f = args.file.name

skip_wraps = 1
strobe_D = get_strobe_events(f, 0x1, skip_wraps=skip_wraps)
strobe_A = get_strobe_events(f, 0x2, skip_wraps=skip_wraps)
delta_D = get_delta_events(f, 0, skip_wraps=skip_wraps)
delta_A = get_delta_events(f, 1, skip_wraps=skip_wraps)

def shifted_deltas(deltas, state, off):
        """ Add an offset to the start times of delta records with the given state """
        ret = np.copy(deltas)
        taken = deltas['state'] == state
        ret['start_t'][taken] += off*delta_clock
        return ret

Dem_Dexc = filter_by_spans(strobe_D, shifted_deltas(delta_D, True, start_exc_offset))
Dem_Aexc = filter_by_spans(strobe_D, shifted_deltas(delta_A, True, start_exc_offset))
Aem_Dexc = filter_by_spans(strobe_A, shifted_deltas(delta_D, True, start_exc_offset))
Aem_Aexc = filter_by_spans(strobe_A, shifted_deltas(delta_A, True, start_exc_offset))

start_t = min(Dem_Dexc['t'][0], Dem_Aexc['t'][0], Aem_Dexc['t'][0], Aem_Aexc['t'][0])
end_t = max(Dem_Dexc['t'][-1], Dem_Aexc['t'][-1], Aem_Dexc['t'][-1], Aem_Aexc['t'][-1])

Dem_Dexc_bins = bin_photons(Dem_Dexc['t'], bin_width*strobe_clock, start_t, end_t)
Dem_Aexc_bins = bin_photons(Dem_Aexc['t'], bin_width*strobe_clock, start_t, end_t)
Aem_Dexc_bins = bin_photons(Aem_Dexc['t'], bin_width*strobe_clock, start_t, end_t)
Aem_Aexc_bins = bin_photons(Aem_Aexc['t'], bin_width*strobe_clock, start_t, end_t)

F_Dem_Dexc = Dem_Dexc_bins['count']
F_Dem_Aexc = Dem_Aexc_bins['count']
F_Aem_Dexc = Aem_Dexc_bins['count']
F_Aem_Aexc = Aem_Aexc_bins['count']

npts = 1000
pl.plot(F_Dem_Dexc[:npts], label='Dem Dexc')
pl.plot(F_Dem_Aexc[:npts], label='Dem Aexc')
pl.plot(F_Aem_Dexc[:npts], label='Aem Dexc')
pl.plot(F_Aem_Aexc[:npts], label='Aem Aexc')
pl.legend()
if args.output is None:
        pl.show()
else:
        pl.savefig('BINS_' + args.output)

#D_F_Dem_Dexc = Dem_Dexc_bins['count']
#F_fret = Aem_Dexc_bins['count'] - Lk - Dir
#A_F_Aem_Aexc = Aem_Aexc_bins['count']

#E_raw_PR = 1. * F_Aem_Dexc / (F_Aem_Dexc + F_Dem_Dexc)
#E_raw_PR = E_raw_PR[np.logical_not(np.isnan(E_raw_PR))]
#S_raw = 1. * (F_Dem_Dexc + F_Aem_Dexc) / (F_Dem_Dexc + F_Aem_Dexc + F_Aem_Aexc)
#S_raw = S_raw[np.logical_not(np.isnan(S_raw))]

def fret_eff(don_bins, acc_bins):
        return 1. * acc_bins / (don_bins + acc_bins)

def plot_fret_eff_hist(ax, donor_bins, acceptor_bins, thresh):
	ctot = donor_bins + acceptor_bins
	t = thresh * mean(ctot)
	take = ctot > t
        td, ta = donor_bins[take], acceptor_bins[take]
        ax.locator_params(nbins=4)
        if len(ta) > 0:
                ax.hist(fret_eff(td, ta), bins=20, histtype='step', range=(0,1))
                ax.set_xlabel('FRET Efficiency')
                ax.set_ylabel('Events')
        ax.text(0.1, 0.75, '$%1.2f bkgrd \/(I > %1.2f \/\mathrm{ppb})$' % (thresh, t), transform=ax.transAxes)

plot_fret_eff_hist(pl.subplot(421), F_Dem_Dexc, F_Aem_Dexc, 0.0)
plot_fret_eff_hist(pl.subplot(422), F_Dem_Dexc, F_Aem_Dexc, 1.0)
plot_fret_eff_hist(pl.subplot(423), F_Dem_Dexc, F_Aem_Dexc, 1.5)
plot_fret_eff_hist(pl.subplot(424), F_Dem_Dexc, F_Aem_Dexc, 2.0)
#plot_fret_eff_hist(pl.subplot(425), F_Dem_Dexc, F_Aem_Dexc, 2.5)
#plot_fret_eff_hist(pl.subplot(426), F_Dem_Dexc, F_Aem_Dexc, 3.0)
#plot_fret_eff_hist(pl.subplot(427), F_Dem_Dexc, F_Aem_Dexc, 5.0)
#plot_fret_eff_hist(pl.subplot(428), F_Dem_Dexc, F_Aem_Dexc, 10.0)

#pl.plot(E_raw_PR[:10000], S_raw[:10000], 'bo')
if args.output is None:
        pl.show()
else:
        pl.savefig('FRET_' + args.output)

