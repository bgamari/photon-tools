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

def load_favia(file):
        data = []
        for l in file:
                d = map(float, l.split(',')[0:5])
                data.append(d)
        data = numpy.array(data)
        times = data[:,0]
        counts = data[:,3]
        var = data[:,4]
        return times, counts, var

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


def fit(times, counts, var):
        # Parameters: [ N, tau_d, a ]
        amp = mean(counts[:5])
        p0 = [ 1/(amp-1), 200e-6, 10 ]

        params, cov_x, infodict, mesg, ier = leastsq(residuals, p0, args=(counts, times, var), full_output=True)
        return params


def fit_single(data):
        times, counts, var = load_favia(data)

        # Subtract out offset
        counts -= 1.0

        # Eliminate early data
        low_time = 1e-6
        var = var[times > low_time]
        counts = counts[times > low_time]
        times = times[times > low_time]

        # Run fit
        params = fit(times, counts, var)
        resid = residuals(params, counts, times, var)
        rel = resid / counts

        # Plot results
        from matplotlib import pyplot as pl
        from mpl_toolkits.axes_grid1 import AxesGrid

        fig = pl.figure(1, (5.,5.))
        grid = AxesGrid(fig, 111, nrows_ncols=(2,1), 
                        axes_pad=0.1, share_all=True, label_mode='L')

        x = linspace(min(times), max(times), 1e6)
        grid[0].semilogx(x, model(params, x), label='Model')
        grid[0].semilogx(times, counts, label='Data', linestyle='None', marker='+')
        grid[0].set_xlabel(r'$\tau$')
        grid[0].set_ylabel(r'$G$')
        grid[0].legend()

        grid[1].errorbar(times, resid, yerr=var)
        grid[1].set_ylabel(r'Fit Residuals')
       
        fig.subplots_adjust(left=0.05, right=0.98)
        pl.draw()
        pl.show()

if __name__ == '__main__':
        fit_single(sys.stdin)


