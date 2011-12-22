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
from matplotlib.gridspec import GridSpec
import matplotlib.text as text
import pickle
import copy
import os.path
import logging
import subprocess
from tempfile import NamedTemporaryFile


def shifted_deltas(deltas, state, off):
	""" Add an offset to the start times of delta records with the given state """
	ret = np.copy(deltas)
	if (off != 0):
		taken = deltas['state'] == state
		ret['start_t'][taken] += off
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
	
	Dem_Dexc = bin_photons(p.Dem_Dexc['t'], bin_width, start_t, end_t)
	Dem_Aexc = bin_photons(p.Dem_Aexc['t'], bin_width, start_t, end_t)
	Aem_Dexc = bin_photons(p.Aem_Dexc['t'], bin_width, start_t, end_t)
	Aem_Aexc = bin_photons(p.Aem_Aexc['t'], bin_width, start_t, end_t)

	return AlexDecomp(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc) 
	
def cut_lead_in(bins):
       	i = 0
       	ddZero = True
       	adZero = True
       	aaZero = True
       	while(ddZero or adZero or aaZero):
       		if(ddZero):
       			if(bins.Dem_Dexc[i]['count'] > 0):
				print "ddZero :", i
       				ddZero = False
       		if(adZero):
       			if(bins.Aem_Dexc[i]['count'] > 0):
       				print "adZero :", i
       				adZero = False
       		if(aaZero):
       			if(bins.Aem_Aexc[i]['count'] > 0):
       				print "aaZero :", i
       				aaZero = False
       		i = i + 1
       	i = i - 1
       	newDD = bins.Dem_Dexc[i:]
       	newDA = bins.Dem_Aexc[i:]
       	newAD = bins.Aem_Dexc[i:]
       	newAA = bins.Aem_Aexc[i:]
		
       	return AlexDecomp(newDD, newDA, newAD, newAA) 



def read_filtered_alex_data(filename, offset):
	alex_photons = None
	
        ddFilename = filename + 'DemDexcFilt'
	daFilename = filename + 'DemAexcFilt'
	adFilename = filename + 'AemDexcFilt'
	aaFilename = filename + 'AemAexcFilt'

	ddFile = None
	daFile = None
	adFile = None
	aaFile = None

	if(os.path.exists(ddFilename)):
		print 'Reading from existing alex photon files'
		
		ddFile = open(ddFilename, 'r')
		daFile = open(daFilename, 'r')
		adFile = open(adFilename, 'r')
		aaFile = open(aaFilename, 'r')
		
		ddPhotons = pickle.load(ddFile)
		daPhotons = pickle.load(daFile)
		adPhotons = pickle.load(adFile)
		aaPhotons = pickle.load(aaFile)

		ddFile.close()
		daFile.close()
		adFile.close()
		aaFile.close()

		alex_photons = AlexDecomp(ddPhotons, daPhotons, adPhotons, aaPhotons)
		
	else:
		print 'Creating new filtered alex photon files'

		skip_wraps = 1
		strobe_D = get_strobe_events(filename, 0x1, skip_wraps=skip_wraps)
		print "strobe_D length ", len(strobe_D)
		strobe_A = get_strobe_events(filename, 0x2, skip_wraps=skip_wraps)
		print "strobe_A length ", len(strobe_A)
		delta_D = get_delta_events(filename, 0, skip_wraps=skip_wraps)[1024:]
		print "delta_D length", len(delta_D)
		delta_A = get_delta_events(filename, 1, skip_wraps=skip_wraps)[1024:]
		print "delta_A length", len(delta_A)
		
		alex_photons = get_alex_photons_ben_bursts(strobe_D, strobe_A, delta_D, delta_A, offset)
		

		ddFile = open(ddFilename, 'w')
		daFile = open(daFilename, 'w')
		adFile = open(adFilename, 'w')
		aaFile = open(aaFilename, 'w')
		
		pickle.dump(alex_photons.Dem_Dexc, ddFile)
		pickle.dump(alex_photons.Dem_Aexc, daFile)
		pickle.dump(alex_photons.Aem_Dexc, adFile)
		pickle.dump(alex_photons.Aem_Aexc, aaFile)
		
		ddFile.close()
		daFile.close()
		adFile.close()
		aaFile.close()
		
	return alex_photons

def plot_alex_burst(ax, tdd, tda, tad, taa, threshold, label, labelHeight, yRange, yTicks, isBottom):
	ax.locator_params(nbins=4)
	e_raw, s_raw = alex_e_s_burst(tdd, tda, tad, taa, threshold)
	#ax.plot(e_raw[:1000000], s_raw[:1000000], 'b+')
	ax.plot(e_raw, s_raw, 'k.')
	ax.set_xlabel('Proximity Ratio')
	ax.set_ylabel('Stoichiometry',  name='Liberation Sans', fontsize=14)
	ax.set_ylim(yRange[0],yRange[1])
	ax.set_yticks(yTicks)

	yTop = ax.get_ylim()[1]
	textY = labelHeight

	ax.text(-0.05, textY, label, name='Liberation Sans', fontsize=14)

	tick_locations = np.arange(6) / 5.0 - 0.1
	ax.set_xticks(tick_locations, minor=False)
	if not isBottom:
		pl.setp(ax.get_xticklabels(), visible=False)
	
	#ax.text(0.1, 0.75, '$%1.2f bkgrd \/(I > %1.2f \/\mathrm{\gamma/bin})$' % (thresh, t), transform=ax.transAxes)
	#pl.ylim(0, 1)
	#pl.xlim(0, 1)
	
	#xtick_locations = (np.arange(6) / 5.0) - 0.1
	#ax.set_xticks(xtick_locations, minor=False)
	
	#ytick_locations = (np.arange(6) / 5.0)
	#ax.set_yticks(ytick_locations, minor=False)
		
	#ax.set_xlim(-0.1, 1.1)
	#ax.set_ylim(-0.1, 1.1)

def alex_e_s_burst(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, threshold, includeZeros=False):
	e = []
	s = []
	for i in range(0, len(Dem_Dexc)):
		a = 1.0 * len(Aem_Dexc[i])
		b = 1.0 * (len(Aem_Dexc[i]) + len(Dem_Dexc[i]))
		c = 1.0 * (len(Dem_Dexc[i]) + len(Aem_Dexc[i]) + len(Aem_Aexc[i]))
		d = 1.0 * (len(Dem_Dexc[i]) + len(Dem_Aexc[i]) + len(Aem_Dexc[i]) + len(Aem_Aexc[i]))
		#if ((d >= 20.0) and (b != 0.0) and (c != 0.0)):
		#if  ((b >= 15.0) and (c >= 15.0) and (d >= 30.0) and (b != 0.0) and (c != 0.0)):
		if  ((b >= threshold) and (b != 0.0) and (c != 0.0)):
			e.append(a / b)
			s.append(b / c)
		else:
			if(includeZeros):
				e.append(-1.0)
				s.append(-1.0)
	return e, s

def fret_eff_burst(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, threshold, includeZeros=False):
	e = []
	for i in range(0, len(Dem_Dexc)):
		a = 1.0 * len(Aem_Dexc[i])
		b = 1.0 * (len(Aem_Dexc[i]) + len(Dem_Dexc[i]))
		c = 1.0 * (len(Dem_Dexc[i]) + len(Aem_Dexc[i]) + len(Aem_Aexc[i]))
		d = 1.0 * (len(Dem_Dexc[i]) + len(Dem_Aexc[i]) + len(Aem_Dexc[i]) + len(Aem_Aexc[i]))
		if  ((b >= threshold) and (b != 0.0) and (c != 0.0)):
			e.append(a / b)
		else:
			if(includeZeros):
				e.append(-1.0)
	return e

def plot_fret_burst(ax, Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, threshold, label, labelHeight, yRange, yTicks, isBottom):
	#ax.locator_params(nbins=4)
	ax.hist(fret_eff_burst(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, threshold), bins=48, histtype='step', range=(-0.1,1.1))
	ax.set_xlabel('Proximity Ratio', name='Liberation Sans', fontsize=14)
	ax.set_ylabel('Events',  name='Liberation Sans', fontsize=14)
	ax.set_ylim(yRange[0], yRange[1])
	ax.set_yticks(yTicks)
	
	yTop = ax.get_ylim()[1]
	textY = labelHeight
	
	ax.text(-0.05, textY, label, fontsize=14)
	
	if not isBottom:
		pl.setp(ax.get_xticklabels(), visible=False)
	
	#ax.text(0.1, 0.75, '$%1.2f bkgrd \/(I > %1.2f \/\mathrm{\gamma/bin})$' % (thresh, t), transform=ax.transAxes)
	#xtick_locations = (np.arange(6) / 5.0) - 0.1
	#ax.set_xticks(xtick_locations, minor=False)
	#ax.set_xlim(-0.1,1.1)

#takes in four streams of filtered photons and a list of bursts
#spits out four streams of lists of bursts
def sort_streams_by_burst(bursts, filtered_photons):
		
	DD = filtered_photons.Dem_Dexc
	DA = filtered_photons.Dem_Aexc
	AD = filtered_photons.Aem_Dexc
	AA = filtered_photons.Aem_Aexc
	

	DDi = 0
	DAi = 0
	ADi = 0
	AAi = 0

	DDBursts = []
	DABursts = []
	ADBursts = []
	AABursts = [] 

	for i in range(0, len(bursts)):
		nextBurst = bursts[i]
		start = nextBurst[0]
		end = nextBurst[1]

		#if i < 50:
		#	print "BURST (", start, ", ", end, ")"
		
		nextDDBurst = []
		while((DDi < len(DD)) and (DD[DDi]['t'] < start)):
		       	#if i < 50:
			#	print "BDDi ", DDi, " DD ", DD[DDi], " start ", start, " end ", end
		       	DDi = DDi+1
		while((DDi < len(DD)) and (DD[DDi]['t'] >= start) and (DD[DDi]['t'] <= end)):
			nextDDBurst.append(DD[DDi])
			#if i < 50:
			#	print "DDi ",DDi, " DD ", DD[DDi]['t'], " start ",start, " end ", end, "nextDDBurst len ", len(DD)
			DDi = DDi+1
		DDBursts.append(nextDDBurst)
		#if i < 50:
		#	print "DDBursts len ", len(DDBursts), "lastBust size ", len(DDBursts[-1]), "nextDDBurst len ", len(nextDDBurst)

		nextDABurst = []
		while((DAi < len(DA)) and (DA[DAi]['t'] < start)):
		       	DAi = DAi+1
		while((DAi < len(DA)) and (DA[DAi]['t'] >= start) and (DA[DAi]['t'] <= end)):
		       	nextDABurst.append(DA[DAi])
		       	#if i < 50:
			#	print "DAi ", DAi, " DA ", DA[DAi]['t'], " start ", start, " end ", end, "nextDABurst len ", len(DA)
			DAi = DAi+1
		DABursts.append(nextDABurst)
		#if i < 50:
		#	print "DABursts len ", len(DABursts), "lastBust size ", len(DABursts[-1]), "nextDABurst len ", len(nextDABurst)

		nextADBurst = []
		while((ADi < len(AD)) and (AD[ADi]['t'] < start)):
		       	ADi = ADi+1
		while((ADi < len(AD)) and (AD[ADi]['t'] >= start) and (AD[ADi]['t'] <= end)):
		       	nextADBurst.append(AD[ADi])
		       	#if i < 50:
			#	print "ADi ",ADi, " AD ", AD[ADi]['t'], " start ",start, " end ", end, "nextADBurst len ", len(AD)
			ADi = ADi+1
		ADBursts.append(nextADBurst)
		#if i < 50:
		#	print "ADBursts len ", len(ADBursts), "lastBust size ", len(ADBursts[-1]), "nextADBurst len ", len(nextADBurst)

		nextAABurst = []
		while((AAi < len(AA)) and (AA[AAi]['t'] < start)):
			AAi = AAi+1
		while((AAi < len(AA)) and (AA[AAi]['t'] >= start) and (AA[AAi]['t'] <= end)):
		       	nextAABurst.append(AA[AAi])
		       	#if i < 50:
			#	print "AAi ",AAi, " AA ", AA[AAi]['t'], " start ",start, " end ", end, "nextAABurst len ", len(AA)
			AAi = AAi+1
		AABursts.append(nextAABurst)
		#if i < 50:
		#	print "AABursts len ", len(AABursts), "lastBust size ", len(AABursts[-1]), "nextAABurst len ", len(nextAABurst)
	
	return AlexDecomp(DDBursts, DABursts, ADBursts, AABursts)
	

def get_alex_photons_ben_bursts(strobe_D, strobe_A, delta_D, delta_A, start_exc_offset=0):
			
	Dem_Dexc = filter_by_spans(strobe_D, shifted_deltas(delta_D, True, start_exc_offset), False)		
	Dem_Aexc = filter_by_spans(strobe_D, shifted_deltas(delta_A, True, start_exc_offset), False)	
	Aem_Dexc = filter_by_spans(strobe_A, shifted_deltas(delta_D, True, start_exc_offset), False)
	Aem_Aexc = filter_by_spans(strobe_A, shifted_deltas(delta_A, True, start_exc_offset), False)
	
	return AlexDecomp(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc) 


def run_ben_burst_alg(filename, bg_rate=1000, burst_rate=4000, burst_length=10):    
	
	bgParam = '--bg-rate='+str(bg_rate)
	brParam = '--burst-rate='+str(burst_rate)
	blParam = '--burst-length='+str(burst_length)
	
	args1 = ['/home/rich/.cabal/bin/bayes-burst-find', bgParam, brParam, blParam, filename]
    
	print "running", args1[0], args1[1], args1[2], args1[3], args1[4]
        
	p1 = subprocess.Popen(args1)
	if p1.wait() != 0:
		print p1.stderr.read()
		raise RuntimeError('Error: exit code %d' % p1.returncode)
    
	return np.loadtxt(f+'.spans', dtype=np.dtype('int64'))

def filter(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh):
	ddBack = mean(Dem_Dexc)
	daBack = mean(Dem_Aexc)
	adBack = mean(Aem_Dexc)
	aaBack = mean(Aem_Aexc)
	
	totalBack = ddBack + adBack
	
	ctot = Dem_Dexc + Aem_Dexc
	t = thresh * totalBack
	take = ctot > t
	tdd, tda, tad, taa = Dem_Dexc[take], Dem_Aexc[take], Aem_Dexc[take], Aem_Aexc[take]
	
	ndd = tdd - ddBack
	nda = tda - daBack
	nad = tad - adBack
	naa = taa - aaBack
	
	return ndd, nda, nad, naa, t


def alex_e_raw(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc):
	return 1.0 * Aem_Dexc / (Aem_Dexc + Dem_Dexc)

def alex_s_raw(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc):
	return 1.0 * (Dem_Dexc + Aem_Dexc) / (Dem_Dexc + Aem_Dexc + Aem_Aexc)

def check_s_raw(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc):
	s_raw = alex_s_raw(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc)
	for i in range(0, len(s_raw)):
		if(s_raw[i] >= 1.0):
			print "raw ", s_raw[i], " DemDexc ", Dem_Dexc[i], " DemAexc ", Dem_Aexc[i], " AemDexc ", Aem_Dexc[i], " AemAexc ", Aem_Aexc[i], " DD + AD ", Dem_Dexc[i] + Aem_Dexc[i], " DD+AD+AA ", Dem_Dexc[i] + Aem_Dexc[i] + Aem_Aexc[i]
	return s_raw
		
		
def plot_alex_raw(ax, Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh, label, labelHeight, yRange, yTicks, isBottom):
	tdd, tda, tad, taa, t = filter(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh)
	ax.locator_params(nbins=4)
	e_raw = alex_e_raw(tdd, tda, tad, taa)
		#s_raw = alex_s_raw(tdd, tda, tad, taa)
	s_raw = check_s_raw(tdd, tda, tad, taa)
	ax.plot(e_raw, s_raw, 'k.')
	ax.set_ylabel('Stoichiometry',  name='Liberation Sans', fontsize=14)
	ax.set_ylim(yRange[0], yRange[1])
	ax.set_yticks(yTicks)
	
	yTop = ax.get_ylim()[1]
	textY = labelHeight
		
	ax.text(-0.05, textY, label, fontsize=14)
	
	tick_locations = np.arange(6) / 5.0 - 0.1
	ax.set_xticks(tick_locations, minor=False)
	if not isBottom:
		pl.setp(ax.get_xticklabels(), visible=False)
	
		

def fret_eff(don_bins, acc_bins):
	return 1. * acc_bins / (don_bins + acc_bins)

def plot_fret_eff_hist(ax, Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh, label, labelHeight, yRange, yTicks, isBottom):
	tdd, tda, tad, taa, t = filter(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh)
	td = tdd
	ta = tad
	#ax.locator_params(nbins=4)
	ax.hist(fret_eff(td, ta), bins=24, histtype='step', range=(-0.1,1.1))
	ax.set_xlabel('Proximity Ratio',  name='Liberation Sans', fontsize=14)
	ax.set_ylabel('Events',  name='Liberation Sans', fontsize=14)
	ax.set_ylim(yRange[0], yRange[1])
	ax.set_yticks(yTicks)
		
	yTop = ax.get_ylim()[1]
	textY = labelHeight
		
	ax.text(-0.05, textY, label, fontsize=14)
		
	if not isBottom:
		pl.setp(ax.get_xticklabels(), visible=False)
	

def get_average(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc):
	ctot = Dem_Dexc + Aem_Dexc
	return mean(ctot)

def plot_bins(ax, bins, color, label, labelHeight, yRange, yTicks, isBottom, isMiddle):
	ax.plot(bins['start_t'], bins['count'], color=color)
		
	isTop = not (isBottom or isMiddle)

	xMin = bins['start_t'][0]
	xMax = bins['start_t'][1000]
		
	ax.set_xlim(xMin, xMax)
	ax.set_ylim(yRange[0], yRange[1])
	ax.set_yticks(yTicks)	

	yTop = ax.get_ylim()[1]
	textY = labelHeight
		
	textX = xMin + 0.05*(xMax - xMin)

	ax.text(textX, textY, label, name='Liberation Sans', fontsize=14)
	if isMiddle:
		#			ax.set_ylabel('Photons per ms',  name='Liberation Sans', fontsize=14)
		ax.set_ylabel('Photons per ms',  name='Liberation Sans', fontsize=14)
	if not isBottom:
		ax.get_xaxis().set_tick_params(labelbottom=False)
	if isBottom:
		ax.set_xlabel('Time (ms)',  name='Liberation Sans', fontsize=14)

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


	photons = get_alex_photons(strobe_D, strobe_A, delta_D, delta_A, delta_clock*start_exc_offset)
	bins = get_alex_bins(photons, strobe_clock*bin_width)
	bins = cut_lead_in(bins)

	jiffyPerMs = 128000

	#zeroBin = (1.0 * bins.Dem_Dexc['start_t'][0]) / jiffyPerMs
	zeroBin = 0

	bins.Dem_Dexc['start_t'] = (1.0 * bins.Dem_Dexc['start_t']) / jiffyPerMs - zeroBin
	bins.Dem_Aexc['start_t'] = (1.0 * bins.Dem_Aexc['start_t']) / jiffyPerMs - zeroBin
	bins.Aem_Dexc['start_t'] = (1.0 * bins.Aem_Dexc['start_t']) / jiffyPerMs - zeroBin
	bins.Aem_Aexc['start_t'] = (1.0 * bins.Aem_Aexc['start_t']) / jiffyPerMs - zeroBin


	F_Dem_Dexc = bins.Dem_Dexc['count']
	F_Dem_Aexc = bins.Dem_Aexc['count']
	F_Aem_Dexc = bins.Aem_Dexc['count']
	F_Aem_Aexc = bins.Aem_Aexc['count']

	#print mean(F_Dem_Dexc), mean(F_Dem_Aexc), mean(F_Aem_Dexc), mean(F_Aem_Aexc)

	#average = get_average(F_Dem_Dexc, F_Dem_Aexc, F_Aem_Dexc, F_Aem_Aexc)

	#threshold = 20.0 / average
	
	pl.figure(figsize=(8,8))

	gs1 = GridSpec(3, 1, width_ratios=[3], height_ratios=[1,1,1])

	ax1 = pl.subplot(gs1[0])
	ax2 = pl.subplot(gs1[1], sharex=ax1)
	ax3 = pl.subplot(gs1[2], sharex=ax1)

	gs1.update(bottom = 0.71, top = 0.98, hspace = 0.001)

	gs2 = GridSpec(2, 1, width_ratios=[3], height_ratios=[4,1])
	
	gs2.update(bottom = 0.08, top = 0.64, hspace = 0.001)
	
	yLabel = 0.4

	plot_bins(ax1, bins.Dem_Dexc, 'k', '(a)', yLabel*16, (0, 16), (4, 8, 12), False, False)
	plot_bins(ax2, bins.Aem_Dexc, 'k', '(b)', yLabel*24, (0, 24), (6, 12, 18), False, True)
	plot_bins(ax3, bins.Aem_Aexc, 'k', '(c)', yLabel*12, (0, 12), (3, 6, 9), True, False)
	
	ax4 = pl.subplot(gs2[0])
	ax5 = pl.subplot(gs2[1], sharex=ax4)
		      
	benSpans = run_ben_burst_alg(f)
	filteredPhotons = read_filtered_alex_data(f, delta_clock*start_exc_offset)
	burstsByStream = sort_streams_by_burst(benSpans, filteredPhotons)
	Bursts_Dem_Dexc = burstsByStream.Dem_Dexc
	Bursts_Dem_Aexc = burstsByStream.Dem_Aexc
	Bursts_Aem_Dexc = burstsByStream.Aem_Dexc
	Bursts_Aem_Aexc = burstsByStream.Aem_Aexc
	
	threshold = 30.0
	
	plot_alex_burst(ax4, Bursts_Dem_Dexc, Bursts_Dem_Aexc, Bursts_Aem_Dexc, Bursts_Aem_Aexc, threshold, '(d)', 0.48, (0.0, 1.05), (0.25, 0.5, 0.75), False)
	plot_fret_burst(ax5, Bursts_Dem_Dexc, Bursts_Dem_Aexc, Bursts_Aem_Dexc, Bursts_Aem_Aexc, threshold, '(e)', 65, (0, 160), (40, 80, 120), True)
	
	#for o in fig.findobj(text.Text):
	#	o.set_fontstyle('italic')

	#for label in ax1.get_xticklabels():
	#	label.set_

        if args.output is None:
		pl.show()
	else:
		pl.savefig(args.output)


