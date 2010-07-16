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
from scipy import *
import scipy.io.array_import
from scipy.optimize import leastsq

def load_data(file):
        data = []
        for l in file:
                d = map(float, l.split(',')[0:5])
                data.append(d)
        return numpy.array(data)

def residuals(p, y, x):
	err = y - peval(x, p)
	return err

def peval(x, p):
	n = p[0]
	tau_d = p[1]
	a = p[2]


	tau_taud = x / tau_d

	b = 1.0 / (1.0 + tau_taud)
	c = 1.0 / (1.0 + (tau_taud / a**2))
	

	return (1.0 / n) * b * sqrt(c)

def derivatives(x, p):
	return [derivativeN(x, p), derivativeTau(x, p), derivativeA(x, p)]

def derivativeN(x, p):
	n = p[0]
	tau_d = p[1]
	a = p[2]

	tau_taud = x / tau_d

	b = 1.0 / (1.0 + tau_taud)
	c = 1.0 / (1.0 + (tau_taud / a**2))
	

	return -(1.0 / n**2) * b * sqrt(c)


def derivativeTau(x, p):
	n = p[0]
	tau_d = p[1]
	a = p[2]

	tau_taud = x / tau_d

	b = 1.0 / (1.0 + tau_taud)
	c = 1.0 / (1.0 + (tau_taud / a**2))
	
	d = b**2 * x / (tau_d**2)
	
	e = sqrt(c**3) * x / (a**2 * tau_d**2)


	return (1.0 / n) * ( (d * sqrt(c)) + (b * e))


def derivativeA(x, p):
	n = p[0]
	tau_d = p[1]
	a = p[2]

	tau_taud = x / tau_d


	b = 1.0 / (1.0 + tau_taud)
	
	c = 1.0 / (1.0 + (tau_taud / a**2))
	
	d = sqrt(c**3) * x / (a**3 * taud)

	return (1.0 / n) *  b * d




data = load_data(sys.stdin)

times = data[:,0]
counts = data[:,3]

# Parameters: [ g, tau_d, a ]
p0 = [0.4, 475.0, 14.0]

params, _ = leastsq(residuals, p0, args=(counts, times), maxfev=2000)

print params
