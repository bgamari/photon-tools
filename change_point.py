#!/usr/bin/python

"""
Implementation of intensity change-point analysis
Ben Gamari, 2011

This follows the work of Yang and Watkins
(J. Phys. Chem. B, Vol 109 (2005), pp.617-628
"""

import numpy as np
from numpy import log

def compute_likelihood(times):
        """ Compute the log-likelihood of a changepoint occurring at each timestep.

            Each element of the returned array L[k] corresponds to the log-likelihood of
            a changepoint occurring at photon times[k].
        """
        N = len(times)
        t = times - times[0]
        T = t[-1]

        k = np.indices(times.shape)
        Vk = t[k] / T
        L = 2*k*log(k/Vk) + 2*(N-k)*log((N-k)/(1-Vk)) - 2*N*log(N)
        return L

def find_tau(times, alpha):
        L = compute_likelihood(times)
        N = len(times)
        np.count(L <= tau) / (N-1) 
        
def find_change_point(times, alpha=0.05):
        for k in range(len(times)):
                pass
                #compute_likelihood(times

def find_change_points(times, alpha=0.05):
        """ Find change points of a data set """
        L = compute_likelihood(times)
        k_max = argmax(L)
        Z = L[k_max]
        tau = compute_tau(1-beta)
        C = (Z-L) <= tau

def test_data():
        times = []
        t = 0
        for i in range(1e5):
                dt = np.random.exponential(1)
                t += dt
                times.append(t)
        for i in range(1e5):
                dt = np.random.exponential(2)
                t += dt
                times.append(t)
        return np.array(times)

if __name__ == '__main__':
        times = test_data()
        L = compute_likelihood(times)[0][1:-1]

        import matplotlib.pyplot as pl
        pl.plot(L)
        pl.show()
