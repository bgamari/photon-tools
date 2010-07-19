#!/usr/bin/python
# vim: set fileencoding=utf-8

# fcs-tools - Tools for FCS data analysis
#
# Copyright © 2010 Ben Gamari
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/ .
#
# Author: Ben Gamari <bgamari@physics.umass.edu>
#

import sys 
import struct
import scipy, numpy
from scipy import sqrt
from scipy.optimize import leastsq
from numpy import min, max, linspace, mean

def load_data(file):
        data = []
        for l in file:
                d = map(float, l.split(',')[0:5])
                data.append(d)
        return numpy.array(data)

def residuals(p, y, x, var):
	err = y - model(p, x)
	return err / var

def model(p, x):
	n = p[0]
	tau_d = p[1]
	a = p[2]

	tau_taud = x / tau_d

	b = 1.0 / (1.0 + tau_taud)
	c = 1.0 / (1.0 + tau_taud * a**-2)
	
	return (1.0 / n) * b * sqrt(c)


data = load_data(sys.stdin)

times = data[:,0]
counts = data[:,3] - 1
var = data[:,4]

# Eliminate early data
low_time = 1e-6
var = var[times > low_time]
counts = counts[times > low_time]
times = times[times > low_time]

# Parameters: [ N, tau_d, a ]
amp = mean(counts[:5])
p0 = [ 1/(amp-1), 200e-6, 10 ]

params, cov_x, infodict, mesg, ier = leastsq(residuals, p0, args=(counts, times, var), full_output=True)
print params, mesg

resid = residuals(params, counts, times, var)
rel = resid / counts

from matplotlib import pyplot as pl
x = linspace(min(times), max(times), 1e6)
pl.semilogx(times, counts, label='Data', linestyle='None', marker='+')
pl.semilogx(x, model(params, x), label='Model')
pl.legend()
pl.show()


