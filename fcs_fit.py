#!/usr/bin/python
# vim: set fileencoding=utf-8

# fcs-tools - Tools for FCS data analysis
#
# Copyright © 2010 Ben Gamari
#
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
	g = p[0]
	tau_d = p[1]
	a = p[2]
	return g/((1 + (x/tau_d)) * (1 + a**-2 * (x/tau_d))**(1/2)) 

def f(x, data):
        res = []
        diff_time, aspect_ratio, nbars = x
        for nbar, (times,counts) in zip(data, nbars):
                tau_taud = times / diff_time
                b = 1.0 + (tau_taud / aspect_ratio**2)
                c = 1.0 / N_bar / (1.0 + tau_taud) / sqrt(b)
                d = abs(counts - c)
                #d /= weight
                res.extend(d)

        return res

data = load_data(sys.stdin)

times = data[:,0]
counts = data[:,1]

# Parameters: [ g, tau_d, a ]
p0 = [5.3, 1.0, 0.3]

params, _ = leastsq(residuals, p0, args=(counts, times))

print params

