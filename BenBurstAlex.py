#!/usr/bin/python

import os.path
import logging
import subprocess
from tempfile import NamedTemporaryFile
import numpy as np
from numpy import log10
from matplotlib import pyplot as pl

import math

from collections import namedtuple
import argparse
from photon_tools.filter_photons import filter_by_spans
from photon_tools.timetag_parse import get_strobe_events, get_delta_events
from photon_tools.bin_photons import bin_photons
from photon_tools.timetag_types import *
from numpy import mean, std, amin, amax, logical_and
import pickle
import copy
from matplotlib.gridspec import GridSpec


def shifted_deltas(deltas, state, off):
	""" Add an offset to the start times of delta records with the given state """
	ret = np.copy(deltas)
	if (off != 0):
		taken = deltas['state'] == state
		ret['start_t'][taken] += off
	return ret

FretDecomp = namedtuple('FretDecomp', 'Dem_Dexc Aem_Dexc')

AlexDecomp = namedtuple('AlexDecomp', 'Dem_Dexc Dem_Aexc Aem_Dexc Aem_Aexc')


def which_span(benSpans, photon):
	for i in range(0,len(benSpans)):
		if((photon >= benSpans[i][0]) and (photon <= benSpans[i][1])):
			return i
	return -1

def check_monotonic_increasing(photon):
	for i in range(1, len(photon)):
		if(photon[i]<photon[i-1]):
			return False
	return True

def count_burst_photons(burstsByStream):
	bursts = len(burstsByStream[0])
	totalPhotons = 0
	
	for i in range(0, bursts):
		totalPhotons += len(burstsByStream[0][i])
		totalPhotons += len(burstsByStream[1][i])
		totalPhotons += len(burstsByStream[2][i])
		totalPhotons += len(burstsByStream[3][i])	
	averagePhotons = 1.0*totalPhotons / bursts		
	return totalPhotons, averagePhotons

def get_alex_values(data):
	return alex_e_s_burst(data[0], data[1], data[2], data[3], True)

def find_values(data, value):
	values = []
	for i in range(0,len(data)):
		if(data[i] == value):
			values.append(i)
	return values

def get_values(values, data):
	valuesOut = []
	for i in values:
		valuesOut.append((len(data[0][i]), len(data[1][i]), len(data[2][i]), len(data[3][i])))
	return valuesOut

def get_data_by_value(data, value):
	alex_values = get_alex_values(data)[0]
	values = find_values(alex_values, value)
	valuesOut = get_values(values, data)
	return valuesOut

def count_burst_sizes(data, n=-1):
	burst_sizes = []
	if(n==-1):
		for i in range(0, len(data[0])):
			burst_sizes.append(len(data[0][i])+len(data[1][i])+len(data[2][i])+len(data[3][i]))
	else:
		for i in range(0, len(data[n])):
			burst_sizes.append(len(data[n][i]))
	return burst_sizes

def find_count(burst_sizes, n):
	count = []
	for i in range(0, len(burst_sizes)):
		if(burst_sizes[i]==n):
			count.append(i)
	return count

def check_for_empty_bursts(burstsByStream):
	bs = count_burst_sizes(burstsByStream)
	return find_count(bs, 0)

def get_smallest_bursts(burstsByStream, burstFloor = 0):
	burst_sizes = count_burst_sizes(burstsByStream)
	minBurst = 10000
	smallest_bursts = []
	for i in range(0, len(burst_sizes)):
		if((burst_sizes[i] < minBurst) and (burst_sizes[i] >= burstFloor)):
			minBurst = burst_sizes[i]
			smallest_bursts = []
			smallest_bursts.append(i)
		elif(burst_sizes[i] == minBurst):
			smallest_bursts.append(i)
	return smallest_bursts

def hist_burst_sizes(burstsByStream):
	burst_sizes = count_burst_sizes(burstsByStream)
	highest = max(burst_sizes)
	hist = []
	for i in range(0, highest+1):
		hist.append(0)
	for i in range(0, len(burst_sizes)):
		hist[burst_sizes[i]] +=1
	return hist

def average_burst(burstsByStream):
	bh = hist_burst_sizes(burstsByStream)
	total = 0
	totalBursts = 0
	for i in range(0, len(bh)):
		total += bh[i]*i
		totalBursts += bh[i]
	return total*1.0/totalBursts

#def skipNPhotons(strobe_D, strobe_A, skip):
#	di = 0
#	ai = 0
#
#	photonCount = 0
#
#	while(photonCount < skip):
#		nextPhoton = min(strobe_D[di][0], strobe_A[ai][0])
#		
#		if(strobe_D[di][0] == nextPhoton):
#			photonCount += 1
#			di += 1
#		elif(strobe_A[ai][0] == nextPhoton):
#			photonCount += 1
#			ai += 1
#		print nextPhoton, ", ", di, ", ", ai
#
#	if(strobe_D[di][0] == nextPhoton):
#		photonCount += 1
#		di += 1
#	if(strobe_A[ai][0] == nextPhoton):
#		photonCount += 1
#		ai += 1
#	
#	nextPhoton = min(strobe_D[di][0], strobe_A[ai][0])
#	
#	return (strobe_D[di, newStrobeA)


def get_alex_photons(strobe_D, strobe_A, delta_D, delta_A, start_exc_offset=0):
		
	Dem_Dexc = filter_by_spans(strobe_D, shifted_deltas(delta_D, True, start_exc_offset))
	Dem_Aexc = filter_by_spans(strobe_D, shifted_deltas(delta_A, True, start_exc_offset))
	Aem_Dexc = filter_by_spans(strobe_A, shifted_deltas(delta_D, True, start_exc_offset))
	Aem_Aexc = filter_by_spans(strobe_A, shifted_deltas(delta_A, True, start_exc_offset))
	return AlexDecomp(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc) 

def filter_strobe_ben_bursts(bursts, strobes):
	filtered_strobes = []
	burstI = 0
	for i in range(0,len(strobes)):
		nextStrobe = strobes[i]
		while(bursts[burstI][1] < nextStrobe[0]):
			burstI += 1
			if(burstI >= len(bursts)):
				return np.array(filtered_strobes, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
		if ((bursts[burstI][0] <= nextStrobe[0]) and (bursts[burstI][1] >= nextStrobe[0])):
			#print "burstI ", burstI, " start ", bursts[burstI][0], "photon ", nextStrobe[0], " end ", bursts[burstI][1]
			filtered_strobes.append(nextStrobe)
		i += 1
	
	return np.array(filtered_strobes, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))

#takes in four streams of filtered photons and a list of bursts
#spits out four streams of lists of bursts
def sort_channel_by_burst(bursts, filtered_photons):
		
	photonI = 0

	outBursts = []
	
	for i in range(0, len(bursts)):
		nextBurst = bursts[i]
		start = nextBurst[0]
		end = nextBurst[1]

		
		nextOutBurst = []
		while((photonI < len(filtered_photons)) and (filtered_photons[photonI]['t'] < start)):
		       	photonI += 1
		while((photonI < len(filtered_photons)) and (filtered_photons[photonI]['t'] >= start) and (filtered_photons[photonI]['t'] <= end)):
			nextOutBurst.append(filtered_photons[photonI])
			photonI += 1
		
		nextOutArray = np.array(nextOutBurst, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))	
		outBursts.append(nextOutArray)
		
	
	return outBursts

	
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

def get_alex_bins(photons, bin_width):
	p = photons
	start_t = min(p.Dem_Dexc['t'][0], p.Dem_Aexc['t'][0], p.Aem_Dexc['t'][0], p.Aem_Aexc['t'][0])
	end_t = max(p.Dem_Dexc['t'][-1], p.Dem_Aexc['t'][-1], p.Aem_Dexc['t'][-1], p.Aem_Aexc['t'][-1])
	
	Dem_Dexc = bin_photons(p.Dem_Dexc['t'], bin_width, start_t, end_t)
	Dem_Aexc = bin_photons(p.Dem_Aexc['t'], bin_width, start_t, end_t)
	Aem_Dexc = bin_photons(p.Aem_Dexc['t'], bin_width, start_t, end_t)
	Aem_Aexc = bin_photons(p.Aem_Aexc['t'], bin_width, start_t, end_t)

	return AlexDecomp(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc) 
	

def merge_streams(DD, DA, AD, AA):
	return np.sort(np.hstack((DD, DA, AD, AA)))

def merge_bins(alexTimeStamps):

	DD = alexTimeStamps.Dem_Dexc
	DA = alexTimeStamps.Dem_Aexc
	AD = alexTimeStamps.Aem_Dexc
	AA = alexTimeStamps.Aem_Aexc
	
	photonCopy = copy.deepcopy(DD)
	for i in range(0, len(DD)):
		photonCopy[i]['count'] += DA[i]['count']
		photonCopy[i]['count'] += AD[i]['count']
		photonCopy[i]['count'] += AA[i]['count']
		
	return photonCopy

def get_alex_bursts(alexTimeStamps, minPhotonsPerBurst = 10, minPhotonsPerWindow = 10, window = 500e-6 * 128e6):	
	
	DD = alexTimeStamps.Dem_Dexc['t']
	DA = alexTimeStamps.Dem_Aexc['t']
	AD = alexTimeStamps.Aem_Dexc['t']
	AA = alexTimeStamps.Aem_Aexc['t']
	
	allPhotons = merge_streams(DD, DA, AD, AA)
	return getBursts(allPhotons, minPhotonsPerBurst, minPhotonsPerWindow, window)

def sortBursts(alexTimeStamps, bursts):
	
	DD = alexTimeStamps.Dem_Dexc
	DA = alexTimeStamps.Dem_Aexc
	AD = alexTimeStamps.Aem_Dexc
	AA = alexTimeStamps.Aem_Aexc
	

	DDi = 0
	DAi = 0
	ADi = 0
	AAi = 0

	DDBursts = []
	DABursts = []
	ADBursts = []
	AABursts = [] 

	DDNoise = []
	DANoise = []
	ADNoise = []
	AANoise = []
	

	for i in range(0, len(bursts)):
		nextBurst = bursts[i]
		start = nextBurst[0]
		end = nextBurst[-1]

		#print "BURST (", start, ", ", end, ")"
		
		nextDDBurst = []
		while((DDi < len(DD)) and (DD[DDi]['t'] < start)):
		       	DDNoise.append(DD[DDi])
		       	#print "BDDi ",DDi, " DD ", DD[DDi], " start ",start, " end ", end
		       	DDi = DDi+1
		while((DDi < len(DD)) and (DD[DDi]['t'] >= start) and (DD[DDi]['t'] <= end)):
			nextDDBurst.append(DD[DDi])
			#print "ADDi ",DDi, " DD ", DD[DDi]['t'], " start ",start, " end ", end
			DDi = DDi+1
		DDBursts.append(nextDDBurst)

		nextDABurst = []
		while((DAi < len(DA)) and (DA[DAi]['t'] < start)):
		       	DANoise.append(DA[DAi])
			DAi = DAi+1
		while((DAi < len(DA)) and (DA[DAi]['t'] >= start) and (DA[DAi]['t'] <= end)):
		       	nextDABurst.append(DA[DAi])
		       	DAi = DAi+1
		DABursts.append(nextDABurst)

		nextADBurst = []
		while((ADi < len(AD)) and (AD[ADi]['t'] < start)):
		       	ADNoise.append(AD[ADi])
		       	ADi = ADi+1
		while((ADi < len(AD)) and (AD[ADi]['t'] >= start) and (AD[ADi]['t'] <= end)):
		       	nextADBurst.append(AD[ADi])
		       	ADi = ADi+1
		ADBursts.append(nextADBurst)

		nextAABurst = []
		while((AAi < len(AA)) and (AA[AAi]['t'] < start)):
			AANoise.append(AA[AAi])
			AAi = AAi+1
		while((AAi < len(AA)) and (AA[AAi]['t'] >= start) and (AA[AAi]['t'] <= end)):
		       	nextAABurst.append(AA[AAi])
		       	AAi = AAi+1
		AABursts.append(nextAABurst)
	
	return (AlexDecomp(DDBursts, DABursts, ADBursts, AABursts), AlexDecomp(DDNoise, DANoise, ADNoise, AANoise))

#this returns lists of photons sorted, instead of keeping everything seperated by bursts
def sortBurstsList(alexTimeStamps, bursts):
	
	DD = alexTimeStamps.Dem_Dexc
	DA = alexTimeStamps.Dem_Aexc
	AD = alexTimeStamps.Aem_Dexc
	AA = alexTimeStamps.Aem_Aexc
	

	DDi = 0
	DAi = 0
	ADi = 0
	AAi = 0

	DDBursts = []
	DABursts = []
	ADBursts = []
	AABursts = []

	DDNoise = []
	DANoise = []
	ADNoise = []
	AANoise = []


	for i in range(0, len(bursts)):
		nextBurst = bursts[i]
		start = nextBurst[0]
		end = nextBurst[-1]

		#print "BURST (", start, ", ", end, ")"
		
		while((DDi < len(DD)) and (DD[DDi]['t'] < start)):
			DDNoise.append(DD[DDi])
			#print "BDDi ",DDi, " DD ", DD[DDi], " start ",start, " end ", end
			DDi = DDi+1
		while((DDi < len(DD)) and (DD[DDi]['t'] >= start) and (DD[DDi]['t'] <= end)):
			DDBursts.append(DD[DDi])
			#print "ADDi ",DDi, " DD ", DD[DDi]['t'], " start ",start, " end ", end
			DDi = DDi+1
		
		while((DAi < len(DA)) and (DA[DAi]['t'] < start)):
		       	DANoise.append(DA[DAi])
		       	DAi = DAi+1
	       	while((DAi < len(DA)) and (DA[DAi]['t'] >= start) and (DA[DAi]['t'] <= end)):
		       	DABursts.append(DA[DAi])
		       	DAi = DAi+1
		
		while((ADi < len(AD)) and (AD[ADi]['t'] < start)):
		       	ADNoise.append(AD[ADi])
		       	ADi = ADi+1
		while((ADi < len(AD)) and (AD[ADi]['t'] >= start) and (AD[ADi]['t'] <= end)):
		       	ADBursts.append(AD[ADi])
		       	ADi = ADi+1
		
		while((AAi < len(AA)) and (AA[AAi]['t'] < start)):
		       	AANoise.append(AA[AAi])
		       	AAi = AAi+1
		while((AAi < len(AA)) and (AA[AAi]['t'] >= start) and (AA[AAi]['t'] <= end)):
		       	AABursts.append(AA[AAi])
		       	AAi = AAi+1
		
	DDBurstsTyped = np.array(DDBursts, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	DABurstsTyped = np.array(DABursts, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	ADBurstsTyped = np.array(ADBursts, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	AABurstsTyped = np.array(AABursts, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))

	DDNoiseTyped = np.array(DDNoise, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	DANoiseTyped = np.array(DANoise, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	ADNoiseTyped = np.array(ADNoise, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	AANoiseTyped = np.array(AANoise, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))


	return (AlexDecomp(DDBurstsTyped, DABurstsTyped, ADBurstsTyped, AABurstsTyped), AlexDecomp(DDNoiseTyped, DANoiseTyped, ADNoiseTyped, AANoiseTyped))

def sortBenBursts(alexTimeStamps, bursts):
	
	DD = alexTimeStamps.Dem_Dexc
	DA = alexTimeStamps.Dem_Aexc
	AD = alexTimeStamps.Aem_Dexc
	AA = alexTimeStamps.Aem_Aexc
	

	DDi = 0
	DAi = 0
	ADi = 0
	AAi = 0

	DDBursts = []
	DABursts = []
	ADBursts = []
	AABursts = [] 

	DDNoise = []
	DANoise = []
	ADNoise = []
	AANoise = []
	

	for i in range(0, len(bursts)):
		nextBurst = bursts[i]
		start = nextBurst[0]
		end = nextBurst[1]

		#print "BURST (", start, ", ", end, ")"
		
		nextDDBurst = []
		while((DDi < len(DD)) and (DD[DDi]['t'] < start)):
		       	DDNoise.append(DD[DDi])
		       	#print "BDDi ",DDi, " DD ", DD[DDi], " start ",start, " end ", end
		       	DDi = DDi+1
		while((DDi < len(DD)) and (DD[DDi]['t'] >= start) and (DD[DDi]['t'] <= end)):
			nextDDBurst.append(DD[DDi])
			#print "ADDi ",DDi, " DD ", DD[DDi]['t'], " start ",start, " end ", end
			DDi = DDi+1
		DDBursts.append(nextDDBurst)

		nextDABurst = []
		while((DAi < len(DA)) and (DA[DAi]['t'] < start)):
		       	DANoise.append(DA[DAi])
			DAi = DAi+1
		while((DAi < len(DA)) and (DA[DAi]['t'] >= start) and (DA[DAi]['t'] <= end)):
		       	nextDABurst.append(DA[DAi])
		       	DAi = DAi+1
		DABursts.append(nextDABurst)

		nextADBurst = []
		while((ADi < len(AD)) and (AD[ADi]['t'] < start)):
		       	ADNoise.append(AD[ADi])
		       	ADi = ADi+1
		while((ADi < len(AD)) and (AD[ADi]['t'] >= start) and (AD[ADi]['t'] <= end)):
		       	nextADBurst.append(AD[ADi])
		       	ADi = ADi+1
		ADBursts.append(nextADBurst)

		nextAABurst = []
		while((AAi < len(AA)) and (AA[AAi]['t'] < start)):
			AANoise.append(AA[AAi])
			AAi = AAi+1
		while((AAi < len(AA)) and (AA[AAi]['t'] >= start) and (AA[AAi]['t'] <= end)):
		       	nextAABurst.append(AA[AAi])
		       	AAi = AAi+1
		AABursts.append(nextAABurst)
	
	return (AlexDecomp(DDBursts, DABursts, ADBursts, AABursts), AlexDecomp(DDNoise, DANoise, ADNoise, AANoise))


#this returns lists of photons sorted, from a list of start and end times
def sortBenBurstsList(alexTimeStamps, bursts):
	
	DD = alexTimeStamps.Dem_Dexc
	DA = alexTimeStamps.Dem_Aexc
	AD = alexTimeStamps.Aem_Dexc
	AA = alexTimeStamps.Aem_Aexc
	

	DDi = 0
	DAi = 0
	ADi = 0
	AAi = 0

	DDBursts = []
	DABursts = []
	ADBursts = []
	AABursts = []

	DDNoise = []
	DANoise = []
	ADNoise = []
	AANoise = []


	for i in range(0, len(bursts)):
		nextBurst = bursts[i]
		start = nextBurst[0]
		end = nextBurst[1]

		#print "BURST (", start, ", ", end, ")"
		
		while((DDi < len(DD)) and (DD[DDi]['t'] < start)):
			DDNoise.append(DD[DDi])
			#print "BDDi ",DDi, " DD ", DD[DDi], " start ",start, " end ", end
			DDi = DDi+1
		while((DDi < len(DD)) and (DD[DDi]['t'] >= start) and (DD[DDi]['t'] <= end)):
			DDBursts.append(DD[DDi])
			#print "ADDi ",DDi, " DD ", DD[DDi]['t'], " start ",start, " end ", end
			DDi = DDi+1
		
		while((DAi < len(DA)) and (DA[DAi]['t'] < start)):
		       	DANoise.append(DA[DAi])
		       	DAi = DAi+1
	       	while((DAi < len(DA)) and (DA[DAi]['t'] >= start) and (DA[DAi]['t'] <= end)):
		       	DABursts.append(DA[DAi])
		       	DAi = DAi+1
		
		while((ADi < len(AD)) and (AD[ADi]['t'] < start)):
		       	ADNoise.append(AD[ADi])
		       	ADi = ADi+1
		while((ADi < len(AD)) and (AD[ADi]['t'] >= start) and (AD[ADi]['t'] <= end)):
		       	ADBursts.append(AD[ADi])
		       	ADi = ADi+1
		
		while((AAi < len(AA)) and (AA[AAi]['t'] < start)):
		       	AANoise.append(AA[AAi])
		       	AAi = AAi+1
		while((AAi < len(AA)) and (AA[AAi]['t'] >= start) and (AA[AAi]['t'] <= end)):
		       	AABursts.append(AA[AAi])
		       	AAi = AAi+1
		
	DDBurstsTyped = np.array(DDBursts, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	DABurstsTyped = np.array(DABursts, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	ADBurstsTyped = np.array(ADBursts, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	AABurstsTyped = np.array(AABursts, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))

	DDNoiseTyped = np.array(DDNoise, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	DANoiseTyped = np.array(DANoise, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	ADNoiseTyped = np.array(ADNoise, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))
	AANoiseTyped = np.array(AANoise, dtype = np.dtype([('t', '<u8'), ('chs', '|u1')]))


	return (AlexDecomp(DDBurstsTyped, DABurstsTyped, ADBurstsTyped, AABurstsTyped), AlexDecomp(DDNoiseTyped, DANoiseTyped, ADNoiseTyped, AANoiseTyped))


def countPhotonsInWindow(photons, i, t):
	min = photons[i] - t/2.0
	max = photons[i] + t/2.0
	index = i-1
	while(photons[index] >= min):
		index -= 1
		if(index < 0):
			break
	count = i - index
	index = i
	while(photons[index] <= max):
		index += 1
		if(index >= len(photons)):
			break
	count += index - i - 1
	return count

def nextBurstStart(photons, i, m, t):
	index = i
	while(countPhotonsInWindow(photons, index, t) < m):
		index += 1
		if(index >= len(photons)):
			index = -1
			break
	return index

def nextBurstEnd(photons, i, m, t):
	index = i
	while(countPhotonsInWindow(photons, index, t) >= m):
		index += 1
		if(index >= len(photons)):
			break
	return index

def getNextBurst(photons, i, l, m, t):
	index = i
	nextBurst = (0, 0)
	while((nextBurst[1] - nextBurst[0]) < l):
		start = nextBurstStart(photons, index, m, t)
		end = nextBurstEnd(photons, start, m, t)
		nextBurst = (start, end)
		if(start < 0):
			nextBurst = (0, 0)
			break
		index = end
	return nextBurst

#see 2006 Weiss "shot noise limited single-molecule....page 3 "definition of bursts"
def getBursts(photons, l, m, t):
	bursts = []
	index = 0
	nextBurst = getNextBurst(photons, 0, l, m, t)
	while((nextBurst[1] - nextBurst[0]) > 0):
		bursts.append(photons[nextBurst[0]:nextBurst[1]])
		index = nextBurst[1]
		nextBurst = getNextBurst(photons, index, l, m, t)
	return bursts

def testGetBursts():
	photons = []
	works = True
	for i in range(0, 1000):
		if (i%10) in (0, 3,4,5,6,7):
			photons.append(i)
	bursts = getBursts(photons, 3, 3, 2.1)
	if(len(bursts) != 100):
		works = False
	count = 0
	for a in bursts:
		count += len(a)
	if(count != 300):
		works = False
	return works

def getPhotonsInBurst(photons, burstMin, burstMax):
	burst = []
	for a in photons:
		if a >= burstMin:
			burst.append(a)
			if a > burstMax:
				break
	return burst


def read_fret_data(filename):
	
	fret_photons = None

	#ddFilename = args.file.name + 'DemDexc'
	#adFilename = args.file.name + 'AemDexc'
	
        ddFilename = f + 'DemDexc'
        adFilename = f + 'AemDexc'

	ddFile = None
	daFile = None

	if(os.path.exists(ddFilename)):
		print 'Reading from existing fret photon files'
		
		ddFile = open(ddFilename, 'r')
		daFile = open(daFilename, 'r')
		
		ddPhotons = pickle.load(ddFile)
		daPhotons = pickle.load(daFile)

		ddFile.close()
		daFile.close()

		fret_photons = FretDecomp(ddPhotons, daPhotons)
		
	else:

		skip_wraps = 1
		strobe_D = get_strobe_events(filename, 0x1, skip_wraps=skip_wraps)[1024:]
		print "Strobe_D length ", len(strobe_D)
		strobe_A = get_strobe_events(filename, 0x2, skip_wraps=skip_wraps)[1024:]
		print "Strobe_A length ", len(strobe_A)
		
		ddFile = open(ddFilename, 'w')
		daFile = open(daFilename, 'w')
		
		pickle.dump(fret_photons.Dem_Dexc, ddFile)
		pickle.dump(fret_photons.Dem_Aexc, daFile)
		
		ddFile.close()
		daFile.close()

		fret_photons = FretDecomp(strobe_D, strobe_A)
		
	return fret_photons
	
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
	
def read_alex_data(filename, offset):
	
	alex_photons = None
	
	#ddFilename = args.file.name + 'DemDexc'
	#daFilename = args.file.name + 'DemAexc'
	#adFilename = args.file.name + 'AemDexc'
	#aaFilename = args.file.name + 'AemAexc'

        ddFilename = filename + 'DemDexc'
	daFilename = filename + 'DemAexc'
	adFilename = filename + 'AemDexc'
	aaFilename = filename + 'AemAexc'

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
		print 'Creating new alex photon files'

		skip_wraps = 1
		strobe_D = get_strobe_events(filename, 0x1, skip_wraps=skip_wraps)[1024:]
		print "strobe_D length ", len(strobe_D)
		strobe_A = get_strobe_events(filename, 0x2, skip_wraps=skip_wraps)[1024:]
		print "strobe_A length ", len(strobe_A)
		delta_D = get_delta_events(filename, 0, skip_wraps=skip_wraps)[1024:]
		print "delta_D length", len(delta_D)
		delta_A = get_delta_events(filename, 1, skip_wraps=skip_wraps)[1024:]
		print "delta_A length", len(delta_A)
		
		alex_photons = get_alex_photons(strobe_D, strobe_A, delta_D, delta_A, offset)
		

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

def shift_alex_data(photons, shiftBy):
	DD = photons.Dem_Dexc
	DA = photons.Dem_Aexc
	AD = photons.Aem_Dexc
	AA = photons.Aem_Aexc

	if(len(DD) + len(DA) + len(AD) + len(AA) <= 1000):
		print "ERROR, Photon Stream Too Short"
		return "ERROR"
	
	
	DDi = 0
	DAi = 0
	ADi = 0
	AAi = 0
	
	photonCount = 0

	while(photonCount < shiftBy):
		nextPhoton = min(DD[DDi]['t'], DA[DAi]['t'], AD[ADi]['t'], AA[AAi]['t'])
		
		if (DD[DDi]['t'] == nextPhoton):
			photonCount += 1
			DDi += 1
		elif (DA[DAi]['t'] == nextPhoton):
			photonCount += 1
			DAi += 1
		elif (AD[ADi]['t'] == nextPhoton):
			photonCount += 1
			ADi += 1
		elif (AA[AAi]['t'] == nextPhoton):
			photonCount += 1
			AAi += 1

		print nextPhoton, ", ", DDi, ", ", DAi, ", ", ADi, ", ", AAi 

       	if (DD[DDi]['t'] == nextPhoton):
		photonCount += 1
		DDi += 1
	if (DA[DAi]['t'] == nextPhoton):
		photonCount += 1
		DDi += 1
	if (AD[ADi]['t'] == nextPhoton):
		photonCount += 1
		DDi += 1
	if (AA[AAi]['t'] == nextPhoton):
		photonCount += 1
		DDi += 1

	nextPhoton = min(DD[DDi]['t'], DA[DAi]['t'], AD[ADi]['t'], AA[AAi]['t'])
	print nextPhoton, ", ", DDi, ", ", DAi, ", ", ADi, ", ", AAi 
	
	return nextPhoton

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

def plot_alex_raw(ax, Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh):
	tdd, tda, tad, taa, t = filter(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh)
	ax.locator_params(nbins=4)
	if len(tdd) > 0:
		e_raw = alex_e_raw(tdd, tda, tad, taa)
		s_raw = alex_s_raw(tdd, tda, tad, taa)
			#ax.plot(e_raw[:1000000], s_raw[:1000000], 'b+')
		ax.plot(e_raw, s_raw, 'k.')
		ax.set_xlabel('Proximity Ratio')
		ax.set_ylabel('Stochiometry')

	#ax.text(0.1, 0.75, '$%1.2f bkgrd \/(I > %1.2f \/\mathrm{\gamma/bin})$' % (thresh, t), transform=ax.transAxes)
	#pl.ylim(0, 1)
	#pl.xlim(0, 1)
	
	xtick_locations = (np.arange(6) / 5.0) - 0.1
	ax.set_xticks(xtick_locations, minor=False)
	
	ytick_locations = (np.arange(6) / 5.0)
	ax.set_yticks(ytick_locations, minor=False)

def alex_e_s_burst(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, threshold, includeZeros=False):
	e = []
	s = []
	for i in range(0, len(Dem_Dexc)):
		a = 1.0 * len(Aem_Dexc[i])
		b = 1.0 * (len(Aem_Dexc[i]) + len(Dem_Dexc[i]))
		c = 1.0 * (len(Dem_Dexc[i]) + len(Aem_Dexc[i]) + len(Aem_Aexc[i]))
		d = 1.0 * (len(Dem_Dexc[i]) + len(Dem_Aexc[i]) + len(Aem_Dexc[i]) + len(Aem_Aexc[i]))
		if  ((b >= threshold) and (b != 0.0) and (c != 0.0)):
			e.append(a / b)
			s.append(b / c)
		else:
			if(includeZeros):
				e.append(-1.0)
				s.append(-1.0)
	return e, s

#def alex_e_burst(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc):
#	e = []
#	for i in range(0, len(Dem_Dexc)):
#		numer = 1.0 * len(Aem_Dexc[i])
#		denom = 1.0 * (len(Aem_Dexc[i]) + len(Dem_Dexc[i]))
#		e.append(numer / denom)
#	return e

#def alex_s_burst(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc):
#	s = []
#	for i in range(0, len(Dem_Dexc)):
#		s.append(1.0 * (len(Dem_Dexc[i]) + len(Aem_Dexc[i])) / (len(Dem_Dexc[i]) + len(Aem_Dexc[i]) + len(Aem_Aexc[i])))
#	return s

def plot_alex_burst(ax, tdd, tda, tad, taa, threshold):
	#ax.locator_params(nbins=4)
	e_raw, s_raw = alex_e_s_burst(tdd, tda, tad, taa, threshold)
	#ax.plot(e_raw[:1000000], s_raw[:1000000], 'b+')
	ax.plot(e_raw, s_raw, 'k.')
	ax.set_xlabel('Proximity Ratio')
	ax.set_ylabel('Stochiometry')

	#ax.text(0.1, 0.75, '$%1.2f bkgrd \/(I > %1.2f \/\mathrm{\gamma/bin})$' % (thresh, t), transform=ax.transAxes)
	#pl.ylim(0, 1)
	#pl.xlim(0, 1)
	
	xtick_locations = (np.arange(6) / 5.0) - 0.1
	ax.set_xticks(xtick_locations, minor=False)
	
	ytick_locations = (np.arange(6) / 5.0)
	ax.set_yticks(ytick_locations, minor=False)
		
	ax.set_xlim(-0.1, 1.1)
	ax.set_ylim(-0.1, 1.1)

def fret_eff(don_bins, acc_bins):
	return 1. * acc_bins / (don_bins + acc_bins)

def plot_fret_eff_hist(ax, Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh):
	tdd, tda, tad, taa, t = filter(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, thresh)
	td = tdd
	ta = tad
	ax.locator_params(nbins=4)
	if len(ta) > 0:
		ax.hist(fret_eff(td, ta), bins=24, histtype='step', range=(-0.1,1.1))
		ax.set_xlabel('Proximity Ratio')
		#ax.set_ylabel('Events')
	ax.text(0.1, 0.75, '$%1.2f bkgrd \/(I > %1.2f \/\mathrm{\gamma/bin})$' % (thresh, t), transform=ax.transAxes)

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

def plot_fret_burst(ax, Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, threshold):
	#ax.locator_params(nbins=4)
	ax.hist(fret_eff_burst(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc, threshold), bins=24, histtype='step', range=(-0.1,1.1))
	#ax.set_xlabel('Proximity Ratio')
	#ax.set_ylabel('Events')
	#ax.text(0.1, 0.75, '$%1.2f bkgrd \/(I > %1.2f \/\mathrm{\gamma/bin})$' % (thresh, t), transform=ax.transAxes)
	
	xtick_locations = (np.arange(6) / 5.0) - 0.1
	ax.set_xticks(xtick_locations, minor=False)
	ax.set_xlim(-0.1,1.1)

def get_average(Dem_Dexc, Dem_Aexc, Aem_Dexc, Aem_Aexc):
	ctot = Dem_Dexc + Aem_Dexc
	return mean(ctot)


def plot_bins(ax, bins, binCount, color):
	ax.plot(bins['start_t'], bins['count'], color=color)
	ax.set_xlim(bins['start_t'][0], bins['start_t'][binCount])
	ax.set_xlabel('Time')
	ax.set_ylabel('Counts')
	
def plot_burst_bins(ax, bins, color):
	return

def plot_burst_hist(ax, bins, color):
	ax.hist(bins['count'], bins=20, log=True, color=color)
	ax.set_xlabel('Burst size (photons)')
	ax.set_ylabel('Events')


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


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	
	parser.add_argument('file', type=argparse.FileType('r'), help='Timetag file')
	parser.add_argument('-s', '--start-offset', type=float, help='Time offset of valid data after delta channel goes high (seconds)', default=1e-6)
	parser.add_argument('-w', '--bin-width', type=float, help='Bin width (seconds)', default=1e-3)
	parser.add_argument('-o', '--output', metavar='FILE', help='Output File Name')

	args = parser.parse_args()

	start_exc_offset = args.start_offset
	bin_width = args.bin_width
	strobe_clock = delta_clock = 128e6
	f = args.file.name
	
        #datapath="/home/rich/data/Sheema/11_08_2011/Alex"
    
        #file1=datapath+"/2011-11-08-run_005.timetag"
    
        #datapath1="/home/rich/data/Alex/09_05_2011"
	
        #file1=datapath1+"/2011-09-05-run_022.timetag"


#############Run Ben's Burst algorithm
	benSpans = run_ben_burst_alg(f)

#############Read in photons normally
	#photons = read_alex_data(f, delta_clock*start_exc_offset)	
        #startingPhoton = shift_alex_data(photons, 1000)

##########Bin them normally
	#bins = get_alex_bins(photons, strobe_clock*bin_width)
	#Bins_Dem_Dexc = bins.Dem_Dexc['count']
	#Bins_Dem_Aexc = bins.Dem_Aexc['count']
	#Bins_Aem_Dexc = bins.Aem_Dexc['count']
	#Bins_Aem_Aexc = bins.Aem_Aexc['count']
	#print "<Bins_Dem_Dexc> ", mean(Bins_Dem_Dexc)," <Bins_Dem_Aexc> ", mean(Bins_Dem_Aexc), " <Bins_Aem_Dexc> ", mean(Bins_Aem_Dexc), " <F_Aem_Aexc> ", mean(Bins_Aem_Aexc)
       	

############Combine the 4 sets of bins and plot them
	#allBins = merge_bins(bins)
	#plot_bins(pl.subplot(613), allBins, 500, 'b')
	
####Align the figure
	gs1 = GridSpec(2, 1, width_ratios=[4], height_ratios=[4,2])
	ax1 = pl.subplot(gs1[0])
	ax2 = pl.subplot(gs1[1])

	#gs1.update(bottom = 0.3, top = 0.95, hspace = 0.001)

	#gs2 = GridSpec(1, 1, width_ratios=[5], height_ratios=[2])
	#ax2 = pl.subplot(gs2[0])

	#gs2.update(bottom = 0.05, top = 0.22, hspace = 0.001)

	threshold = 30.0

####Plot the Alex using Ben's Burst detection
	#get the photons from the two files, without removing the empty spacing
	filteredPhotons = read_filtered_alex_data(f, delta_clock*start_exc_offset)	
	burstsByStream = sort_streams_by_burst(benSpans, filteredPhotons)
	totalBurstPhotons, averageBurstPhotons = count_burst_photons(burstsByStream)
	print "Total Burst Photons", totalBurstPhotons, "Average Photons Per Burst", averageBurstPhotons
	Bursts_Dem_Dexc = burstsByStream.Dem_Dexc
	Bursts_Dem_Aexc = burstsByStream.Dem_Aexc
	Bursts_Aem_Dexc = burstsByStream.Aem_Dexc
	Bursts_Aem_Aexc = burstsByStream.Aem_Aexc
	plot_alex_burst(ax1, Bursts_Dem_Dexc, Bursts_Dem_Aexc, Bursts_Aem_Dexc, Bursts_Aem_Aexc, threshold)
	
####Plot the Proximity Ratio using Ben's Burst detection
	plot_fret_burst(ax2, Bursts_Dem_Dexc, Bursts_Dem_Aexc, Bursts_Aem_Dexc, Bursts_Aem_Aexc, threshold)

###########Combine the 4 sets of bins and plot them with the 50 microsecond gaps
	#spreadBins = get_alex_bins(filteredPhotons, strobe_clock*bin_width*2)
	#Spread_Bins_Dem_Dexc = spreadBins.Dem_Dexc['count']
	#Spread_Bins_Dem_Aexc = spreadBins.Dem_Aexc['count']
	#Spread_Bins_Aem_Dexc = spreadBins.Aem_Dexc['count']
	#Spread_Bins_Aem_Aexc = spreadBins.Aem_Aexc['count']
	#allSpreadBinsSpaced = merge_bins(spreadBins)	
	#plot_bins(pl.subplot(612), allSpreadBinsSpaced, 500, 'b')
	
###########Combine the 4 sets of ben's burst filtered bins and plot the used and not used
	#usedPhotons, noisePhotons = sortBurstsList(filteredPhotons, benSpans)

	#usedBins = get_alex_bins(usedPhotons, strobe_clock*bin_width)
		
	#allUsedBinsSpaced = merge_bins(usedBins)
	
	#noiseBins = get_alex_bins(noisePhotons, strobe_clock*bin_width)

	#allNoiseBinsSpaced = merge_bins(noiseBins)

	#plot_bins(pl.subplot(613), allUsedBinsSpaced, 500, 'b')

	#plot_bins(pl.subplot(614), allNoiseBinsSpaced, 500, 'b')

	
############Calculate the average number of photons per bin, and some reasonable potential thresholds to use.
	
	#average = get_average(Bins_Dem_Dexc, Bins_Dem_Aexc, Bins_Aem_Dexc, Bins_Aem_Aexc)
	
	#thresholdOne = 0.0
	#thresholdTwo = 5.0 / average
	#thresholdThree = 2.0 * thresholdTwo
	#thresholdFour = 3.0 * thresholdTwo
	#thresholdFive = 4.0 * thresholdTwo

####Plot the Alex using Rich's burst detection, paramater set one
        
    #burstWindows = get_alex_bursts(photons, 10, 10, 500e-6*128e6)
    
    #bursts = sortBurstsList(photons, burstWindows)[0]
    
    
    
    #burstBins = get_alex_bins(bursts, strobe_clock*bin_width)
    
    #allBurstBins = merge_bins(burstBins)
    
    
    #plot_bins(pl.subplot(434), allBurstBins, 500, 'r')
    
    #burstsSorted = sortBursts(photons, burstWindows)[0]
    
    #Bursts_Dem_Dexc = burstsSorted.Dem_Dexc
    #Bursts_Dem_Aexc = burstsSorted.Dem_Aexc
    #Bursts_Aem_Dexc = burstsSorted.Aem_Dexc
    #Bursts_Aem_Aexc = burstsSorted.Aem_Aexc
    
    #plot_alex_burst(pl.subplot(435), Bursts_Dem_Dexc, Bursts_Dem_Aexc, Bursts_Aem_Dexc, Bursts_Aem_Aexc)
    
  
####Plot the Alex using Rich's burst detection paramater set two

    #burstWindows = get_alex_bursts(photons, 6, 6, 500e-6*128e6)
    
    #bursts = sortBurstsList(photons, burstWindows)[0]
    
    
    #burstBins = get_alex_bins(bursts, strobe_clock*bin_width)
    
    #allBurstBins = merge_bins(burstBins)
    
    
    #plot_bins(pl.subplot(437), allBurstBins, 500, 'r')
    
    #burstsSorted = sortBursts(photons, burstWindows)[0]
    
    #Bursts_Dem_Dexc = burstsSorted.Dem_Dexc
    #Bursts_Dem_Aexc = burstsSorted.Dem_Aexc
    #Bursts_Aem_Dexc = burstsSorted.Aem_Dexc
    #Bursts_Aem_Aexc = burstsSorted.Aem_Aexc
    
    #plot_alex_burst(pl.subplot(438), Bursts_Dem_Dexc, Bursts_Dem_Aexc, Bursts_Aem_Dexc, Bursts_Aem_Aexc)
    

  
####Plot the Alex using Rich's burst detection paramater set 3
   
    
    #burstWindows = get_alex_bursts(photons, 2, 2, 500e-6*128e6)
    
    #bursts = sortBurstsList(photons, burstWindows)[0]
    
    
    #burstBins = get_alex_bins(bursts, strobe_clock*bin_width)
    
    #allBurstBins = merge_bins(burstBins)
    
    
    #plot_bins(pl.subplot(4,3,10), allBurstBins, 500, 'r')
    
    
    #burstsSorted = sortBursts(photons, burstWindows)[0]
    
    #Bursts_Dem_Dexc = burstsSorted.Dem_Dexc
    #Bursts_Dem_Aexc = burstsSorted.Dem_Aexc
    #Bursts_Aem_Dexc = burstsSorted.Aem_Dexc
    #Bursts_Aem_Aexc = burstsSorted.Aem_Aexc
    
    #plot_alex_burst(pl.subplot(4,3,11), Bursts_Dem_Dexc, Bursts_Dem_Aexc, Bursts_Aem_Dexc, Bursts_Aem_Aexc)
    

#####Output everything to either the screen or a file

	if args.output is None:
		pl.show()
	else:
		pl.savefig(args.output)


