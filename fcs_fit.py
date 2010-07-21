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
import collections
DataSet = collections.namedtuple('DataSet', 'name times counts var')

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

def model(p, x):
	tau_d = p[0]
	a = p[1]
	n = p[2]

	tau_taud = x / tau_d

	b = 1.0 / (1.0 + tau_taud)
	c = 1.0 / (1.0 + tau_taud * a**-2)
	
	return (1.0 / n) * b * sqrt(c)

def fitfunc(p, data):
        res = []
        for d, n in zip(data, p[2:]):
                p1 = list(p[0:2]) + [n]
                err = d.counts - model(p1, d.times)
                res.extend(err / d.var)
	return numpy.array(res)

def fit(data):
        est_n = lambda d: 1/(mean(d.counts[:5]) - 1)
        p0 = [ 200e-6, 10 ] + map(est_n, data)
        params, cov_x, infodict, mesg, ier = leastsq(fitfunc, p0, args=(data), full_output=True)
        return params

data = []
for f in sys.argv[1:]:
        times, counts, var = load_favia(open(f))

        # Subtract out offset
        counts -= 1.0

        # Eliminate early data
        low_time = 1e-6
        var = var[times > low_time]
        counts = counts[times > low_time]
        times = times[times > low_time]

        data.append(DataSet(f, times, counts, var))

# Run fit
params = fit(data)

# Plot results
from matplotlib import pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable

ax = pl.subplot(111)
x = linspace(min(times), max(times), 1e6)

for d,n in zip(data, params[2:]):
        p1 = list(params[0:2]) + [n]
        ax.semilogx(x, model(p1, x), label='%s (fit, N=%f)' % (d.name, n))
        ax.semilogx(d.times, d.counts, label=d.name, linestyle='None', marker='+')

divider = make_axes_locatable(ax)
for d,n in zip(data, params[2:]):
        p1 = list(params[0:2]) + [n]
        resid = model(p1, d.times) - d.counts

        ax_resid = divider.append_axes("top", 1.4, pad=0.0, sharex=ax)
        ax_resid.axhline(0, color='black')
        ax_resid.errorbar(times, resid, yerr=var, linestyle='None', marker='x')
        pl.setp(ax_resid.get_xticklabels(), visible=False)
        ax_resid.set_ylabel(r'Fit Residuals')

text = [
        r'$\tau_d = \mathrm{%1.3e}$' % params[0],
        r'$a = \mathrm{%1.3e}$' % params[1],
]
pl.figtext(0.70, 0.40, "\n".join(text))

ax.set_xlabel(r'$\tau$')
ax.set_ylabel(r'$G$')
ax.legend()

ax.autoscale_view(tight=True, scalex=True)
pl.draw()
pl.show()


#if __name__ == '__main__':
#        print fit_single(sys.stdin)


