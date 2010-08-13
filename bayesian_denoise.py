#!/usr/bin/python

import numpy as np
import pywt
from math import log

def cutoff_denoise(bins, cutoff, level=1):
        a,*d = pywt.wavedec(bins, 'haar', level=level)
        d[d > thresh] = 0
        return pywt.waverec([a] + d, 'haar', level=level)

def denoise(bins, cutoff, level=1):
        n = len(d1)
        sigma = bins**0.5 # Poisson shot noise
        tau = sigma * (2*log(len(bins)))**0.5
        return cutoff_denoises(bins, tau)

