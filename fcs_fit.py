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
from numpy import min, max, linspace

def load_data(file):
        data = []
        for l in file:
                d = map(float, l.split(',')[0:5])
                data.append(d)
        return numpy.array(data)

def residuals(p, y, x):
	err = y - model(p, x)
	return err

def model(p, x):
	n = p[0]
	tau_d = p[1]
	a = p[2]
        offset = p[3]

	tau_taud = x / tau_d

	b = 1.0 / (1.0 + tau_taud)
	c = 1.0 / (1.0 + tau_taud * a**-2)
	
	return (1.0 / n) * b * sqrt(c) + offset


data = load_data(sys.stdin)

times = data[:,0]
counts = data[:,3]

# Parameters: [ g, tau_d, a, offset ]
p0 = [1,1,1,0.1]

params, cov_x, infodict, mesg, ier = leastsq(residuals, p0, args=(counts, times), full_output=True)
print params, mesg, ier

resid = residuals(params, counts, times)
rel = resid / counts

from matplotlib import pyplot as pl
x = linspace(min(times), max(times))
pl.plot(times, counts, label='Data')
pl.plot(x, model(params, x), label='Model')
pl.legend()
pl.show()

