#!/usr/bin/python

""" 
Implementation of FRET denoising using wavelet decomposition
Ben Gamari, 2010

This follows the work of Taylor, et al.
(Biophysical Journal, Vol 98, January 2010, pp. 164-173)
"""

import numpy as np
from numpy import abs, mean, std, sign
import pywt
from math import log
import copy
import random

def hard_threshold(bins, cutoff, level=1):
        coeffs = pywt.wavedec(bins, 'haar', level=level)
        np.place(coeffs[1:], abs(coeffs[1:]) > cutoff, 0)
        return pywt.waverec([a,d], 'haar')

def soft_threshold(bins, tau, level=1, plot=False):
        coeffs = pywt.wavedec(bins, 'haar', level=level)
        old_coeffs = copy.deepcopy(coeffs) if plot else None
        coeffs[1:] = [sign(a) * np.maximum(0, abs(a)-tau) for a in coeffs[1:] ]
        filtered = pywt.waverec(coeffs, 'haar')

        if plot:
                from matplotlib import pyplot as pl
                fig = pl.figure(figsize=(15,8))
                fig.suptitle("Level %d Wavelet Decomposition" % level)
                fig.subplots_adjust(hspace=0.3, wspace=0.3)

                ax = fig.add_subplot(2,1,1)
                ax.plot(bins, label='Raw')
                ax.plot(filtered, label='Filtered')
                ax.legend()

                ax = fig.add_subplot(2,level,1+level)
                ax.set_title("Approximation coefficients")
                ax.plot(coeffs[0])
                for l in range(1,level):
                        ax = fig.add_subplot(2, level, 1+l+level)
                        ax.plot(old_coeffs[l])
                        ax.axhline(-tau, color='g'); ax.axhline(+tau, color='g')
                        ax.set_title("Level %d detail coefficients" % l)

                fig.savefig('wavelet-analysis.png')

        return filtered


def denoise(bins, level=1, plot=False):
        sigma = np.mean(bins)**0.5 # Poisson shot noise
        tau = sigma * (2*log(len(bins)))**0.5
        return soft_threshold(bins, tau, level, plot)

def test_data(transitions=100):
        states = [ 163, 220, 230, 100, 40, 45 ]

        clean, noisy = [], []
        for i in range(transitions):
                length = random.randint(10, 100)
                state = random.choice(states)
                clean.append([state] * length)
                noisy.append(np.random.normal(state, 17, length))
                
        return np.hstack(clean), np.hstack(noisy)

if __name__ == "__main__":
        level = 4

        noisy = np.random.normal(163, 17, 500)
        denoised = denoise(noisy, level)
        print 'raw', mean(noisy), std(noisy)
        print 'denoised', mean(denoised), std(denoised)

        clean, noisy = test_data()
        denoised = denoise(noisy, level)
        #err = (noisy - denoised)**2

        from matplotlib import pyplot as pl
        pl.clf()
        pl.plot(noisy, label='noisy', alpha=0.3)
        pl.plot(denoised, label='denoised')
        pl.plot(clean, label='clean')
        pl.legend()
        pl.show()

