#!/usr/bin/python

""" 
Implementation of Bayesian blink removal algorithm.
Ben Gamari, 2010

This follows the work of Taylor, et al.
(Biophysical Journal, Vol 98, January 2010, pp. 164-173)

"""

import logging
from math import pi, exp
import random
import numpy as np
from numpy import mean, array

logging.basicConfig(level=logging.DEBUG)

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
                px = fret_bins.A - bgA
                pd = fret_bins.D - bgD + px
                ct_param = np.mean(px / pd)
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
                _max = max(self.fret_bins.D)
                _min = min(self.fret_bins.D)
                center = (_max - _min)/2 + _min
                blinks = self.fret_bins.D < center
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
                        if True:
                                rng = 100000
                                from matplotlib import pyplot as pl
                                from mpl_toolkits.axes_grid1 import make_axes_locatable
                                pl.clf()
                                ax = pl.subplot(111)
                                ax.plot(self.fret_bins.A[:rng], label='Acceptor')
                                ax.plot(self.fret_bins.D[:rng], label='Donor')
                                ax.legend()
                                divider = make_axes_locatable(ax)
                                cax = divider.append_axes('bottom', size='100%', pad=0.2)
                                cax.plot(B[:rng], label='P(B|Na)')
                                cax.plot(NB[:rng], label='P(NB|Na)')
                                cax.legend()
                                cax = divider.append_axes('bottom', size='100%', pad=0.4)
                                cax.plot(bayes[:rng], label='Bayes factor')
                                cax.set_ylim(0, bayes_thresh*4)
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
def test_data(transitions=1e2, bg_flux=10, flux=(220, 10), fret_eff=0.40, ct_prob=0.005):
        """ Produce fake FRET data. """
        fret_length = int(round(transitions*0.75))
        ct_length = int(round(transitions*0.25))
        bg_length = transitions * 5e2
        d_blink_prob = 0.2
        a_blink_prob = 0.2

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

        # FRET region
        fret_a_bins, fret_d_bins = [], []
        for i in range(fret_length):
                if random.random() < d_blink_prob:
                        # Donor blink
                        length = random.randint(0, 5*100)
                        fret_d_bins.append(np.zeros(length))
                        fret_a_bins.append(np.zeros(length))
                elif random.random() < a_blink_prob:
                        # Acceptor blink
                        length = random.randint(0, 5*100)
                        fret_d_bins.append(np.random.normal(flux[0], flux[1], (length)).round())
                        fret_a_bins.append(np.zeros(length))
                else:
                        # No blinking
                        length = int(random.normalvariate(1e4, 1e3))
                        fret_d_bins.append(np.random.normal((1-fret_eff)*flux[0], flux[1], (length)).round())
                        fret_a_bins.append(np.random.normal(fret_eff*flux[0], flux[1], (length)).round())

        fret_bins = np.rec.fromarrays([np.hstack(fret_a_bins), np.hstack(fret_d_bins)], names='A,D')
        noisify_bins(fret_bins)

        # Crosstalk region (acceptor died)
        ct_a_bins, ct_d_bins = [], []
        for i in range(ct_length):
                if random.random() < d_blink_prob:
                        # Donor blink
                        length = random.randint(0, 5*100)
                        ct_d_bins.append(np.zeros(length))
                        ct_a_bins.append(np.zeros(length))
                else:
                        # No blinking
                        length = round(random.normalvariate(1e4, 1e3))
                        ct_d_bins.append(np.random.normal(flux[0], flux[1], (length)))
                        ct_a_bins.append(np.zeros(length))

        ct_bins = np.rec.fromarrays([np.hstack(ct_a_bins), np.hstack(ct_d_bins)], names='A,D')
        noisify_bins(ct_bins)

        # Background region
        bg_d_bins = np.random.normal(bg_flux, bg_flux, bg_length).round()
        bg_a_bins = np.random.normal(bg_flux, bg_flux, bg_length).round()
        bg_bins = np.rec.fromarrays([bg_a_bins, bg_d_bins], names='A,D')
        np.place(bg_bins.A, bg_bins.A<0, 0)
        np.place(bg_bins.D, bg_bins.D<0, 0)

        return fret_bins, ct_bins, bg_bins

if __name__ == '__main__':
        bayes_thresh = 1e6

        from matplotlib import pyplot as pl
        fret,ct,bg = test_data(transitions=1e2)
        bins = np.hstack((fret,ct,bg))

        pl.clf()
        pl.plot(bins['A'], label='Acceptor')
        pl.plot(bins['D'], label='Donor')
        pl.legend()
        pl.savefig('initial.png')

        traj = fret_trajectory.from_bins(fret, ct, bg)
        blinks,_,_ = traj.find_A_blinks(bayes_thresh)
        deblinked = fret[~blinks]

        pl.clf()
        pl.plot(fret.D, label='Donor')
        pl.plot(fret.A, label='Acceptor')
        ymin, ymax = pl.ylim()
        pl.vlines(np.nonzero(blinks)[0], ymin, ymax, alpha=0.3, color='r', label='Blinks')
        pl.legend()
        pl.savefig('fret.png')
        
