#!/usr/bin/python

""" 
Implementation of Bayesian blink removal algorithm.
Ben Gamari, 2010

This follows the work of Taylor, et al.
(Biophysical Journal, Vol 98, January 2010, pp. 164-173)

"""

from math import pi, exp
import random
import numpy as np
from numpy import mean, array

dtype = np.dtype({'names': ['A','D'], 'formats': [int,int]})
class fret_trajectory:
        """ Represents a FRET trajectory. This consists of two-channel photon
            bin series for three regions: the FRET region, the cross-talk region,
            and the background region. These bin series are passed as numpy
            arrays consisting of records with two fields, A (the acceptor
            channel) and D (the donor channel) """

        def __init__(self, fret_bins, ct_bins, bg_bins):
                self.fret_bins = fret_bins
                self.fret_mean = mean(fret_bins)
                self.ct_bins = ct_bins
                self.ct_mean = mean(ct_bins)
                self.bg_bins = bg_bins
                self.bg_mean = mean(bg_bins)
                self.bins = fret_bins + ct_bins + bg_bins

                px = fret_bins['A'] - self.bg_mean['A']
                pd = fret_bins['D'] - self.bg_mean['D'] + px
                self.ct_param = px / pd

        def prior_NB_prob(self, Na):
                """ Calculate P(Na|NB) """
                mu = self.bg_mean['A']
                return (2*pi*mu)**-0.5 * exp(-(Na-mu)**2 / 2 / mu)

        def prior_B_prob(self, Na):
                """ Calculate P(Na|B) """
                mu = self.bg_mean['A'] + self.ct_param * self.fret_mean['D']
                return (2*pi*mu)**-0.5 * exp(-(Na-mu)**2 / 2 / mu)

        def post_B_prob(self, PB, PNB, Na):
                """ Calculate P(B|Na) """
                y = self.prior_B_prob(Na)*PB + prior_NB_prob(Na)*PNB
                return self.prior_B_prob(Na) * PB / y

        def post_NB_prob(self, PB, PNB, Na):
                """ Calculate P(NB|Na) """
                y = self.prior_B_prob(Na)*PB + prior_NB_prob(Na)*PNB
                return self.prior_NB_prob(Na) * PNB / y

        def find_blinks(self, bayes_thresh=2.0):
                old_PB = None
                PB = 0.001
                PNB = 1 - PB
                while True:
                        B = post_B_prob(PB, PNB, self.bins)
                        NB = post_NB_prob(PB, PNB, self.bins)
                        bayes = B / NB
                        blinks = bayes > bayes_thresh
                        nblinks = len(np.nonzero(blinks))
                        fb = nblinks / (len(self.bins) - nblinks)
                        old_PB = PB
                        PB, PNB = fb, 1-fb
                        if abs(PB-old_PB) / old_PB < 0.05:
                                return PB, PNB, blinks

        def remove_blinks(self, bayes_thresh=2.0):
                _, _, blinks = self.find_blinks(bayes_thresh)
                return self.bins[not blinks]

def test_data(transitions=1e3, bg_flux = (10, 5), flux = (220, 10), fret_eff = 0.15, ct_prob=0.01):
        a_bins = []
        d_bins = []

        # FRET region
        for i in range(transitions*0.75):
                state = random.choice(['obs', 'd_blink', 'a_blink'])
                if state == 'obs':
                        length = int(random.normalvariate(1e4, 1e3))
                        d_bins.extend(np.random.normal((1-fret_eff)*flux[0], flux[1], (length)))
                        a_bins.extend(np.random.normal(fret_eff*flux[0], flux[1], (length)))
                elif state == 'd_blink':
                        length = random.randint(0, 5*100)
                        d_bins.extend([0]*length)
                        a_bins.extend([0]*length)
                elif state == 'a_blink':
                        length = random.randint(0, 5*100)
                        d_bins.extend(np.random.normal(fret_eff*flux[0], flux[1], (length)))
                        a_bins.extend([0]*length)

        # Crosstalk region (acceptor died)
        for i in range(transitions*0.10):
                state = random.choice(['obs', 'd_blink'])
                if state == ['obs']:
                        length = int(random.normalvariate(1e4, 1e3))
                        d_bins.extend(np.random.normal(flux[0], flux[1], (length)))
                        a_bins.extend([0]*length)
                elif state == 'd_blink':
                        length = random.randint(0, 5*100)
                        d_bins.extend([0]*length)
                        a_bins.extend([0]*length)

        # Background region
        d_bins.extend([0]*int(0.15*transitions))
        a_bins.extend([0]*int(0.15*transitions))
                
        # Simulate crosstalk
        crosses = np.random.randint(-5, +5, len(d_bins))
        d_bins -= crosses
        a_bins += crosses

        # Add background 
        a_bins += np.random.normal(bg_flux[0], bg_flux[1])
        d_bins += np.random.normal(bg_flux[0], bg_flux[1])
        return array([a_bins, d_bins], dtype=dtype)

if __name__ == '__main__':
        import matplotlib.pyplot as pl
        bg_flux = (10, 5)
        data = test_data(bg_flux=bg_flux)
        #denoised = denoise(data)
        pl.plot(data)
        pl.plot(denoised)
        pl.savefig('hi.png')
        
