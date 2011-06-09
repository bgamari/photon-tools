#!/usr/bin/python

""" 
Implementation of Bayesian blink removal algorithm.
Ben Gamari, 2010

This follows the work of Taylor, et al.
(Biophysical Journal, Vol 98, January 2010, pp. 164-173)
"""

import logging
from math import pi, exp, log
import random
import numpy as np
from numpy import mean, array
from numpy.lib.recfunctions import stack_arrays 

logging.basicConfig(level=logging.DEBUG)
plot_iterations = True
plot_len = 4000

dtype = np.dtype({'names': ['A','D'], 'formats': [np.uint32,np.uint32]})
class fret_trajectory(object):
        """ Represents a FRET trajectory. This consists of two-channel photon
            bin series for three regions: the FRET region, the cross-talk region,
            and the background region. These bin series are passed as numpy
            arrays consisting of records with two fields, A (the acceptor
            channel) and D (the donor channel) """

        def __init__(self, fret_bins, bgA, bgD, ct_param):
                self.fret_bins = fret_bins
                self.fretD_mean = mean(fret_bins.D)
                self.fretA_mean = mean(fret_bins.A)
                self.bgA = bgA
                self.bgD = bgD
                self.ct_param = ct_param
                self.ct_photons = self.ct_param / (1-self.ct_param) * (self.fret_bins.D - self.bgD) # nx
                self.nD = self.fret_bins.D - self.bgD + self.ct_photons
                self.nA = self.fret_bins.A - self.bgA - self.ct_photons

        @classmethod
        def from_bins(klass, fret_bins, ct_bins, bg_bins):
                """ Create a FRET trajectory from bins from FRET, cross-talk,
                    and background regions """
                bgA = np.mean(bg_bins.A)
                bgD = np.mean(bg_bins.D)
                px = ct_bins.A - bgA
                pd = ct_bins.D - bgD + px
                ct_param = np.mean(px / pd)
                logging.debug("Cross-talk parameter=%f" % ct_param)
                return klass(fret_bins, bgA, bgD, ct_param)

        def prior_NB_prob(self, Na):
                """ Calculate P(Na|NB) """
                mu = self.fretA_mean
                return (2*pi*mu)**-0.5 * np.exp(-(Na-mu)**2 / 2 / mu)

        def prior_B_prob(self, Na):
                """ Calculate P(Na|B) """
                mu = self.bgA + self.ct_param * self.ct_photons
                return (2*pi*mu)**-0.5 * np.exp(-(Na-mu)**2 / 2 / mu)

        def post_B_prob(self, PB, PNB, Na):
                """ Calculate P(B|Na) """
                y = self.prior_B_prob(Na)*PB + self.prior_NB_prob(Na)*PNB
                return self.prior_B_prob(Na) * PB / y

        def post_NB_prob(self, PB, PNB, Na):
                """ Calculate P(NB|Na) """
                y = self.prior_B_prob(Na)*PB + self.prior_NB_prob(Na)*PNB
                return self.prior_NB_prob(Na) * PNB / y

        def find_D_blinks(self):
                """ Simple heuristic to identify donor blinks """
                Dmax = max(self.fret_bins.D)
                Dmin = min(self.fret_bins.D)
                Dcenter = (Dmax - Dmin)/2 + Dmin
                Amax = max(self.fret_bins.A)
                Amin = min(self.fret_bins.A)
                Acenter = (Amax - Amin)/2 + Amin
                blinks = (self.fret_bins.D < Dcenter) & (self.fret_bins.A < Acenter)
                return blinks

        def find_A_blinks(self, bayes_thresh=2.0):
                """ Identify acceptor channel blinking with Bayesian inference """
                old_PB = None
                PB = 0.001
                PNB = 1 - PB
                i=0
                while True:
                        B = self.post_B_prob(PB, PNB, self.fret_bins.A)
                        NB = self.post_NB_prob(PB, PNB, self.fret_bins.A)
                        bayes = B / NB
                        if plot_iterations:
                                from matplotlib import pyplot as pl
                                from mpl_toolkits.axes_grid1 import make_axes_locatable
                                pl.clf()
                                ax = pl.subplot(111)
                                ax.plot(self.fret_bins.A[:plot_len], label='Acceptor')
                                ax.plot(self.fret_bins.D[:plot_len], label='Donor')
                                ax.legend()
                                divider = make_axes_locatable(ax)
                                cax = divider.append_axes('bottom', size='100%', pad=0.2)
                                cax.plot(B[:plot_len], label='P(B|Na)')
                                cax.plot(NB[:plot_len], label='P(NB|Na)')
                                cax.legend()
                                cax = divider.append_axes('bottom', size='100%', pad=0.4)
                                x = range(len(bayes[:plot_len]))
                                cax.fill_between(x, bayes[:plot_len], color='b')
                                cax.axhline(bayes_thresh, color='r')
                                cax.set_ylim(0, 8*bayes_thresh)
                                pl.savefig('iter%d.pdf' % i)

                        blinks = bayes > bayes_thresh
                        n_blinks = len(np.nonzero(blinks)[0])
                        fb = 1.0 * n_blinks / self.fret_bins.shape[0]
                        logging.debug('Blink fraction %f' % fb)
                        if abs(PB-fb) / PB < 0.05:
                                return blinks, PB, PNB
                        PB, PNB = fb, 1-fb
                        i += 1

        def find_blinks(self, bayes_thresh=2.0):
                blinksA,_,_ = self.find_A_blinks(bayes_thresh)
                #blinksD = self.find_D_blinks()
                return blinksA #| blinksD

        def remove_blinks(self, bayes_thresh=2.0):
                return self.fret_bins[~self.find_blinks()]


###
### Test Code
###

def kinetic_mc(states, steps):
        """ Generate simulated FRET trajectory by Kinetic Monte Carlo method
        
            States are specified as a dictionary in the form of,
              states = { 'state_name': (rate, event_func), ... }
         
            For every step, a state is chosen and the corresponding event_func
            is called. The state trajectory is returned as a list of state names.
        """
        state_traj = []
        total_rate = sum(rate for rate,f in states.values())
        for i in range(steps):
                s = random.random()*total_rate
                for state, (rate,f) in states.items():
                        if rate > s:
                                state_traj.append(state)
                                dt = int(round(-log(random.random()) / total_rate))
                                f(dt)
                                break
                        else:
                                s -= rate

        return state_traj

def test_data(transitions=1e4, bg_flux=10, flux=220, fret_eff1=0.40, fret_eff2=0.7, ct_prob=0.005):
        """ Produce fake FRET data.
        
            We generate data for each of the three regions of the experiment
            (FRET, crosstalk, background). The FRET region data is generated
            with two FRET states of the given efficiencies.
        """
        fret_length = int(round(transitions*0.5))
        ct_length = int(round(transitions*0.5))
        bg_length = 10000

        bins = []
        add_bins = lambda a,d: bins.append(np.rec.fromarrays((a,d), names='A,D'))

        # FRET region
        #       State   : (rate, event_func)
        states = {
                'Ablink': (2e-3, lambda l: add_bins(np.zeros(l, dtype='i8'),
                                                    np.random.poisson(flux, l).round())),
                'Dblink': (2e-3, lambda l: add_bins(np.zeros(l, dtype='i8'),
                                                    np.zeros(l, dtype='i8'))),
                'fret1' : (6e-3, lambda l: add_bins(np.random.poisson((1-fret_eff1)*flux, l).round(),
                                                    np.random.poisson(fret_eff1*flux, l).round())),
                'fret2' : (8e-3, lambda l: add_bins(np.random.poisson((1-fret_eff2)*flux, l).round(),
                                                    np.random.poisson(fret_eff2*flux, l).round())),
        }
        state_traj = kinetic_mc(states, fret_length)
        fret_bins = stack_arrays(bins, asrecarray=True)

        def noisify_bins(bins):
                # Cross-talk
                crosses = (bins.D * np.random.normal(ct_prob, ct_prob/10, bins.shape)).round()
                bins.D -= crosses
                bins.A += crosses
                # Background
                bins.A += np.random.normal(bg_flux, bg_flux, bins.shape).round()
                bins.D += np.random.normal(bg_flux, bg_flux, bins.shape).round()
                np.place(bins.A, bins.A<0, 0)
                np.place(bins.D, bins.D<0, 0)

        noisify_bins(fret_bins)

        # Crosstalk region (acceptor died)
        bins = []
        states = {
                'Dblink': (2e-3, lambda l: add_bins(np.zeros(l, dtype='i8'),
                                                    np.zeros(l, dtype='i8'))),
                'obs'   : (8e-3, lambda l: add_bins(np.zeros(l, dtype='i8'), 
                                                    np.random.poisson(flux, l).round())),
        }
        state_traj = kinetic_mc(states, fret_length)
        ct_bins = stack_arrays(bins, asrecarray=True)
        noisify_bins(ct_bins)

        # Background region
        bg_d_bins = np.random.poisson(bg_flux, bg_length).round()
        bg_a_bins = np.random.poisson(bg_flux, bg_length).round()
        bg_bins = np.rec.fromarrays([bg_a_bins, bg_d_bins], names='A,D')
        np.place(bg_bins.A, bg_bins.A<0, 0)
        np.place(bg_bins.D, bg_bins.D<0, 0)

        return fret_bins, ct_bins, bg_bins

if __name__ == '__main__':
        bayes_thresh = 10
        transitions = 1e4

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('-t', '--test', action='store_true', help='Use generated test data')
        parser.add_argument('-a', '--acceptor', type=argparse.FileType('r'), help='Acceptor bin data')
        parser.add_argument('-d', '--donor', type=argparse.FileType('r'), help='Donor bin data')
        parser.add_argument('-f', '--fret', metavar='START:STOP', help='FRET region')
        parser.add_argument('-c', '--crosstalk', metavar='START:STOP', help='Cross-talk region')
        parser.add_argument('-b', '--background', metavar='START:STOP', help='Background region')
        args = parser.parse_args()
        
        fret,ct,bg = None,None,None
        if args.test:
                logging.info("Generating test data")
                from matplotlib import pyplot as pl
                fret,ct,bg = test_data(transitions=transitions)
                bins = stack_arrays((fret,ct,bg))
                logging.info("Generated %d time steps" % len(bins))
        else:
                #da = np.fromfile(args.acceptor, dtype='u8,u2')[:,1]
                #dd = np.fromfile(args.donor, dtype='u8,u2')[:,1]
                da = np.fromfile(args.acceptor, dtype='<u2')
                dd = np.fromfile(args.donor, dtype='<u2')
                bins = np.rec.fromarrays([da, dd], names='A,D')
                def parse_range(s):
                        start,stop = s.split(':')
                        return slice(int(start), int(stop))
                fret = bins[parse_range(args.fret)]
                ct = bins[parse_range(args.crosstalk)]
                bg = bins[parse_range(args.background)]

        from matplotlib import pyplot as pl
        pl.clf()
        pl.plot(bins['D'], label='Donor')
        pl.plot(bins['A'], label='Acceptor')
        if args.test:
                pl.axvspan(0, len(fret), color='b', alpha=0.1)
                pl.axvspan(len(fret), len(ct), color='g', alpha=0.1)
                pl.axvspan(len(fret)+len(ct), len(fret)+len(ct)+len(bg), color='r', alpha=0.1)
        else:
                def plot_range(rng, color):
                        r = parse_range(rng)
                        pl.axvspan(r.start, r.stop, color=color, alpha=0.1)
                plot_range(args.fret, 'b')
                plot_range(args.crosstalk, 'g')
                plot_range(args.background, 'r')
        pl.legend()
        pl.savefig('all.pdf')

        pl.clf()
        pl.plot(fret.D[:plot_len], label='Donor')
        pl.plot(fret.A[:plot_len], label='Acceptor')
        pl.legend()
        pl.savefig('initial.pdf')

        logging.info("Finding acceptor blinks")
        traj = fret_trajectory.from_bins(fret, ct, bg)
        blinks,_,_ = traj.find_A_blinks(bayes_thresh)
        deblinked = fret[~blinks]

        pl.clf()
        pl.plot(fret.D[:plot_len], label='Donor')
        pl.plot(fret.A[:plot_len], label='Acceptor')
        ymin, ymax = pl.ylim()
        b = np.nonzero(blinks)[0]
        pl.vlines(b[b<plot_len], ymin, ymax, alpha=0.1, color='r', label='Acceptor blinks')
        pl.legend()
        pl.savefig('fret.pdf')
        
        pl.clf()
        pl.plot(deblinked['D'][:plot_len], label='Donor')
        pl.plot(deblinked['A'][:plot_len], label='Acceptor')
        pl.legend()
        pl.savefig('deblinked.pdf')

