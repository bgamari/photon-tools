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
class fret_trajectory:
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
                Dmax = max(self.fret_bins.D)
                Dmin = min(self.fret_bins.D)
                Dcenter = (Dmax - Dmin)/2 + Dmin
                Amax = max(self.fret_bins.A)
                Amin = min(self.fret_bins.A)
                Acenter = (Amax - Amin)/2 + Amin
                blinks = (self.fret_bins.D < Dcenter) & (self.fret_bins.A < Acenter)
                return blinks

        def find_A_blinks(self, bayes_thresh=2.0):
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
                                pl.savefig('iter%d.png' % i)

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
                blinksD = self.find_D_blinks()
                return blinksA | blinksD

        def remove_blinks(self, bayes_thresh=2.0):
                return self.fret_bins[~self.find_blinks()]


###
### Test Code
###
def test_data(transitions=1e3, bg_flux=10, flux=(220, 10), fret_eff1=0.40, fret_eff2=0.7, ct_prob=0.005):
        """ Produce fake FRET data. """
        fret_length = int(round(transitions*0.75))
        ct_length = int(round(transitions*0.25))
        bg_length = 10000

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

        def kinetic_mc(states, steps):
                state_traj = []
                total_rate = sum(rate for rate,f in states.values())
                for i in range(steps):
                        s = random.random()*total_rate
                        for state, (rate,f) in states.items():
                                if rate > s:
                                        state_traj.append(state)
                                        dt = int(round(-log(random.random()) / total_rate))
                                        #print state, dt
                                        f(dt)
                                        break
                                else:
                                        s -= rate


                return state_traj


        bins = []
        add_bins = lambda a,d: bins.append(np.rec.fromarrays((a,d), names='A,D'))

        # FRET region
        #       State   : (rate, event_func)
        states = {
                'Ablink': (2e-3, lambda l: add_bins(np.zeros(l),
                                                    np.random.normal(flux[0], flux[1], l).round())),
                'Dblink': (2e-3, lambda l: add_bins(np.zeros(l),
                                                    np.zeros(l))),
                'fret1' : (8e-3, lambda l: add_bins(np.random.normal((1-fret_eff1)*flux[0], flux[1], l).round(),
                                                    np.random.normal(fret_eff1*flux[0], flux[1], l).round())),
                'fret2' : (8e-3, lambda l: add_bins(np.random.normal((1-fret_eff2)*flux[0], flux[1], l).round(),
                                                    np.random.normal(fret_eff2*flux[0], flux[1], l).round())),
        }
        state_traj = kinetic_mc(states, fret_length)
        fret_bins = stack_arrays(bins, asrecarray=True)
        noisify_bins(fret_bins)

        # Crosstalk region (acceptor died)
        bins = []
        states = {
                'Dblink': (2e-3, lambda l: add_bins(np.zeros(l),
                                                    np.zeros(l))),
                'obs'   : (8e-3, lambda l: add_bins(np.zeros(l), 
                                                    np.random.normal(flux[0], flux[1], l).round())),
        }
        state_traj = kinetic_mc(states, fret_length)
        ct_bins = stack_arrays(bins, asrecarray=True)
        noisify_bins(ct_bins)

        # Background region
        bg_d_bins = np.random.normal(bg_flux, bg_flux, bg_length).round()
        bg_a_bins = np.random.normal(bg_flux, bg_flux, bg_length).round()
        bg_bins = np.rec.fromarrays([bg_a_bins, bg_d_bins], names='A,D')
        np.place(bg_bins.A, bg_bins.A<0, 0)
        np.place(bg_bins.D, bg_bins.D<0, 0)

        return fret_bins, ct_bins, bg_bins

if __name__ == '__main__':
        bayes_thresh = 10

        from matplotlib import pyplot as pl
        fret,ct,bg = test_data(transitions=1e3)
        bins = stack_arrays((fret,ct,bg))

        pl.clf()
        pl.plot(fret.D[:plot_len], label='Donor')
        pl.plot(fret.A[:plot_len], label='Acceptor')
        pl.legend()
        pl.savefig('initial.png')

        traj = fret_trajectory.from_bins(fret, ct, bg)
        blinksA,_,_ = traj.find_A_blinks(bayes_thresh)
        blinksD = traj.find_D_blinks()
        deblinked = fret[~blinksA & ~blinksD]
        #deblinked = fret[~blinksA]

        pl.clf()
        pl.plot(fret.D[:plot_len], label='Donor')
        pl.plot(fret.A[:plot_len], label='Acceptor')
        ymin, ymax = pl.ylim()
        b = np.nonzero(blinksD)[0]
        pl.vlines(b[b<plot_len], ymin, ymax, alpha=0.1, color='y', label='Donor blinks')
        b = np.nonzero(blinksA)[0]
        pl.vlines(b[b<plot_len], ymin, ymax, alpha=0.1, color='r', label='Acceptor blinks')
        pl.legend()
        pl.savefig('fret.png')
        
        pl.clf()
        pl.plot(deblinked['D'][:plot_len], label='Donor')
        pl.plot(deblinked['A'][:plot_len], label='Acceptor')
        pl.legend()
        pl.savefig('deblinked.png')

