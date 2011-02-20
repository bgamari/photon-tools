#!/usr/bin/python

"""
Implementation of intensity change-point analysis
Ben Gamari, 2011

This follows the work of Yang and Watkins
(J. Phys. Chem. B, Vol 109 (2005), pp.617-628
"""

import numpy as np

def test_data(states, npts):
        """ Generate simulated trajectory of given length with a given set of states """
        d = np.empty(npts)
        i = 0
        state = 0
        while i < npts:
                n,rate = states[state]
                dt = np.random.poisson(n)
                if i+dt >= npts:
                        d.resize(i)
                        return d
                d[i:i+dt] = np.random.poisson(rate, dt)
                i += dt

def compute_likelihood(times):
        """ Compute the log-likelihood of a changepoint occurring at each timestep.

            Each element of the returned array L[k] corresponds to the log-likelihood of
            a changepoint occurring at photon times[k].
        """
        N = len(data)
        t = times - times[0]
        T = t[-1]

        k = np.indices(photons.shape)
        Vk = t[k] / T
        L = 2*k*log(k/Vk) + 2*(N-k)*log((N-k)/(1-Vk)) - 2*N*log(N)
        return L

def find_tau(times, alpha):
        L = compute_likelihood(times)
        N = len(times)
        np.count(L <= tau) / (N-1) - 
        
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


if __name__ == '__main__':
        #d = test_data(states, 1e5)

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
        times = np.array(times)

        L = compute_likelihood(times)[0][1:-1]

        import matplotlib.pyplot as pl
        pl.plot(L)
        pl.show()
