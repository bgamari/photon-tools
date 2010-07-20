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

def fitfunc(p, y, x, var):
	err = y - model(p, x)
	return err / var

def model(p, x):
	tau_d = p[0]
	a = p[1]
	n = p[2]

	tau_taud = x / tau_d

	b = 1.0 / (1.0 + tau_taud)
	c = 1.0 / (1.0 + tau_taud * a**-2)
	
	return (1.0 / n) * b * sqrt(c)


def fit(times, counts, var):
        # Parameters: [ tau_d, a, N ]
        amp = mean(counts[:5])
        p0 = [ 200e-6, 10, 1/(amp-1) ]

        params, cov_x, infodict, mesg, ier = leastsq(fitfunc, p0, args=(counts, times, var), full_output=True)
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
        resid = counts - model(params, times)

        # Plot results
        from matplotlib import pyplot as pl
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        pl.title(data.name)
        ax = pl.subplot(111)

        x = linspace(min(times), max(times), 1e6)
        ax.semilogx(x, model(params, x), label='Fit')
        ax.semilogx(times, counts, label=data.name, linestyle='None', marker='+')
        ax.set_xlabel(r'$\tau$')
        ax.set_ylabel(r'$G$')
        ax.legend()

        text = [
                r'$\tau_d = \mathrm{%1.3e}$' % params[0],
                r'$a = \mathrm{%1.3e}$' % params[1],
                r'$N = \mathrm{%1.3e}$' % params[2],
        ]
        pl.figtext(0.70, 0.40, "\n".join(text))

        divider = make_axes_locatable(ax)
        ax_resid = divider.append_axes("top", 1.4, pad=0.0, sharex=ax)
        ax_resid.axhline(0, color='black')
        ax_resid.errorbar(times, resid, yerr=var, linestyle='None', marker='x')
        pl.setp(ax_resid.get_xticklabels(), visible=False)
        ax_resid.set_ylabel(r'Fit Residuals')
       
        ax.autoscale_view(tight=True, scalex=True)
        pl.draw()
        pl.show()

        return params

if __name__ == '__main__':
        print fit_single(sys.stdin)


