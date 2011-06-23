#!/usr/bin/python

from collections import namedtuple
import argparse
import numpy as np
from photon_tools.filter_photons import filter_by_spans
from photon_tools.timetag_parse import get_strobe_events, get_delta_events
from photon_tools.bin_photons import bin_photons
from photon_tools.timetag_types import *
from matplotlib import pyplot as pl
from numpy import mean, std, amin, amax, logical_and

def shifted_deltas(deltas, state, off):
	""" Add an offset to the start times of delta records with the given state """
	ret = np.copy(deltas)
	if (off != 0):
		taken = deltas['state'] == state
		ret['start_t'][taken] += off*delta_clock
	return ret

AlexDecomp = namedtuple('AlexDecomp', 'Dem_Dexc Dem_Aexc Aem_Dexc Aem_Aexc')

def get_alex_photons(strobe_D, strobe_A, delta_D, delta_A, start_exc_offset=0):
	Dem_Dexc = filter_by_spans(strobe_D, shifted_deltas(delta_D, True, start_exc_offset))
	Dem_Aexc = filter_by_spans(strobe_D, shifted_deltas(delta_A, True, start_exc_offset))
	Aem_Dexc = filter_by_spans(strobe_A, shifted_deltas(delta_D, True, start_exc_offset))
	Aem_Aexc = filter_by_spans(strobe_A, shifted_deltas(delta_A, True, start_exc_offset))
	return AlexDecomp(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc) 

def get_alex_bins(photons, bin_width):
	p = photons
	start_t = min(p.Dem_Dexc['t'][0], p.Dem_Aexc['t'][0], p.Aem_Dexc['t'][0], p.Aem_Aexc['t'][0])
	end_t = max(p.Dem_Dexc['t'][-1], p.Dem_Aexc['t'][-1], p.Aem_Dexc['t'][-1], p.Aem_Aexc['t'][-1])
	
	Dem_Dexc= bin_photons(p.Dem_Dexc['t'], bin_width, start_t, end_t)
	Dem_Aexc= bin_photons(p.Dem_Aexc['t'], bin_width, start_t, end_t)
	Aem_Dexc= bin_photons(p.Aem_Dexc['t'], bin_width, start_t, end_t)
	Aem_Aexc= bin_photons(p.Aem_Aexc['t'], bin_width, start_t, end_t)

	return AlexDecomp(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc) 
	

if __name__ == '__main__':
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
	strobe_D = get_strobe_events(f, 0x1, skip_wraps=skip_wraps)[1024:]
	strobe_A = get_strobe_events(f, 0x2, skip_wraps=skip_wraps)[1024:]
	delta_D = get_delta_events(f, 0, skip_wraps=skip_wraps)[1024:]
	delta_A = get_delta_events(f, 1, skip_wraps=skip_wraps)[1024:]

	photons = get_alex_photons(strobe_D, strobe_A, delta_D, delta_A, start_exc_offset)
	bins = get_alex_bins(photons, strobe_clock*bin_width)
	F_Dem_Dexc = bins.Dem_Dexc['count']
	F_Dem_Aexc = bins.Dem_Aexc['count']
	F_Aem_Dexc = bins.Aem_Dexc['count']
	F_Aem_Aexc = bins.Aem_Aexc['count']

	npts = 10000
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


	def alex_e_raw(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc):
		return 1.0 * Aem_Dexc / (Aem_Dexc + Dem_Dexc)

	def alex_s_raw(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc):
		return 1.0 * (Dem_Dexc + Aem_Dexc) / (Dem_Dexc + Aem_Dexc + Aem_Aexc)

	def plot_alex_raw(ax, Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh):
		ctot = Dem_Dexc + Dem_Aexc + Aem_Dexc + Aem_Aexc
		t = thresh * mean(ctot)
		take = ctot > t
		tdd, tda, tad, taa = Dem_Dexc[take], Dem_Aexc[take], Aem_Dexc[take], Aem_Aexc[take]
		ax.locator_params(nbins=4)
		if len(tdd) > 0:
			e_raw = alex_e_raw(tdd, tda, tad, taa)
			s_raw = alex_s_raw(tdd, tda, tad, taa)
			ax.plot(e_raw[:10000], s_raw[:10000], 'b+')

		ax.text(0.1, 0.75, '$%1.2f bkgrd \/(I > %1.2f \/\mathrm{ ppb})$' % (thresh, t), transform=ax.transAxes)
		pl.ylim(0, 1)
		pl.xlim(0, 1)
		tick_locations = np.arange(10) / 5.
		ax.set_xticks(tick_locations, minor=False)
		ax.set_yticks(tick_locations, minor=False)
		

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
		ax.text(0.1, 0.75, '$%1.2f bkgrd \/(I > %1.2f \/\mathrm{ ppb})$' % (thresh, t), transform=ax.transAxes)

	plot_fret_eff_hist(pl.subplot(241), F_Dem_Dexc, F_Aem_Dexc, 0.0)
	plot_fret_eff_hist(pl.subplot(242), F_Dem_Dexc, F_Aem_Dexc, 1.0)
	plot_fret_eff_hist(pl.subplot(243), F_Dem_Dexc, F_Aem_Dexc, 1.5)
	plot_fret_eff_hist(pl.subplot(244), F_Dem_Dexc, F_Aem_Dexc, 2.0)
	plot_alex_raw(pl.subplot(245), F_Dem_Dexc, F_Dem_Aexc, F_Aem_Dexc, F_Aem_Aexc, 0.0)
	plot_alex_raw(pl.subplot(246), F_Dem_Dexc, F_Dem_Aexc, F_Aem_Dexc, F_Aem_Aexc, 1.0)
	plot_alex_raw(pl.subplot(247), F_Dem_Dexc, F_Dem_Aexc, F_Aem_Dexc, F_Aem_Aexc, 1.5)
	plot_alex_raw(pl.subplot(248), F_Dem_Dexc, F_Dem_Aexc, F_Aem_Dexc, F_Aem_Aexc, 2.0)


	#pl.plot(E_raw_PR[:10000], S_raw[:10000], 'bo')
	if args.output is None:
		pl.show()
	else:
		pl.savefig('FRET_' + args.output)


