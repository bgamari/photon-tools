#!/usr/bin/python

from photon_tools.timetag_parse import get_strobe_events
from photon_tools.bin_photons import bin_photons
import argparse
from numpy import mean, std, amin, amax, logical_and
from matplotlib import pyplot as pl

parser = argparse.ArgumentParser()
parser.add_argument('file', type=file, help='Time trace file')
parser.add_argument('-s', '--bin-size', metavar='TIME', help='Length of histogram bins (seconds)', default=1e-3)
parser.add_argument('-S', '--skip', metavar='N', type=int, help='Skip first N records', default=0)
parser.add_argument('-c', '--clockrate', metavar='HZ', type=float, help='Clockrate of timetagger (will read from .params file by default)', default=None)
parser.add_argument('-A', '--acceptor', metavar='N', required=True, type=int, help='Acceptor channel')
parser.add_argument('-D', '--donor', metavar='N', required=True, type=int, help='Donor channel')
parser.add_argument('-o', '--output', metavar='FILE', help='Output File Name')

args = parser.parse_args()

da = get_strobe_events(args.file.name, 1<<(args.acceptor-1))[args.skip:]
dd = get_strobe_events(args.file.name, 1<<(args.donor-1))[args.skip:]
clockrate = args.clockrate # TODO: Read from params

ba = bin_photons(da['t'], args.bin_size*clockrate)
bd = bin_photons(dd['t'], args.bin_size*clockrate)

# Make sure data are aligned
start_t = max(amin(ba['start_t']), amin(bd['start_t']))
end_t = min(amax(ba['start_t']), amax(bd['start_t']))
ba = ba[logical_and(ba['start_t'] >= start_t, ba['start_t'] < end_t)]
bd = bd[logical_and(bd['start_t'] >= start_t, bd['start_t'] < end_t)]

ctot = ba['count'] + bd['count']
cavg = mean(ctot)
cavg = 0.7

def threshold(thresh):
        take = ctot > thresh*cavg
        return (ba[take], bd[take])

def fret_eff(acc_bins, don_bins):
        return 1. * acc_bins['count'] / (don_bins['count']+acc_bins['count'])

pl.figure()
pl.subplots_adjust(hspace=0.4, left=0.1)

def plot_bins(ax, bins):
        ax.plot(bins['start_t'], bins['count'])
        ax.set_xlim(bins['start_t'][0], bins['start_t'][1000])
        ax.set_xlabel('Time')
        ax.set_ylabel('Counts')

def plot_burst_hist(ax, bins):
        ax.hist(bins['count'], bins=20, log=True)
        ax.set_xlabel('Burst size (photons)')
        ax.set_ylabel('Events')

plot_bins(pl.subplot(421), bd)
plot_burst_hist(pl.subplot(422), ba)

plot_bins(pl.subplot(423), bd)
plot_burst_hist(pl.subplot(424), bd)

def plot_fret_eff_hist(ax, thresh):
        ta,td = threshold(thresh)
        if len(ta) > 0:
                ax.hist(fret_eff(ta, td), bins=20, histtype='step', range=(0,1), log=True)
                ax.set_xlabel('FRET Efficiency')
                ax.set_ylabel('Events')
		ax.text(0.1, 0.75, '$%1.2f \sigma$' % thresh, transform=ax.transAxes)

plot_fret_eff_hist(pl.subplot(425), 1.0)
plot_fret_eff_hist(pl.subplot(426), 1.5)
plot_fret_eff_hist(pl.subplot(427), 2.0)
plot_fret_eff_hist(pl.subplot(428), 4.0)

if args.output is None:
	pl.show()
else:
	pl.savefig(args.output)
