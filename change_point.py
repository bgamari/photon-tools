#!/usr/bin/python

"""
Implementation of intensity change-point analysis
Ben Gamari, 2011

This follows the work of Yang and Watkins
(J. Phys. Chem. B, Vol 109 (2005), pp.617-628
"""

import numpy as np
from numpy import log, sum, sqrt

def tau95(n):
        if 1 < n <= 100:
                return 1.*n/(3.0123 + 0.4835*log(n) - 0.00957*log(n**2) - 0.001488*log(n**3))
        elif 100 < n:
                return 1.*n/(3.0806 + 0.4894*log(n) - 0.02086*log(n**2))
        else:
                raise RuntimeError('Implement me')

def compute_likelihood(times):
        """ Compute the log-likelihood of a changepoint occurring at each timestep.

            Each element of the returned array L[k] corresponds to the log-likelihood of
            a changepoint occurring at photon times[k].
        """
        N = len(times)
        T = times[-1]

        k = np.indices(times.shape)[0]
        Vk = times / T
        L = 2*k*log(k/Vk) + 2*(N-k)*log((N-k)/(1-Vk)) - 2*N*log(N)
        return L

def find_tau(times, alpha):
        L = compute_likelihood(times)
        N = len(times)
        np.count(L <= tau) / (N-1) 
        
def find_change_points(times, alpha=0.05):
        """ Find change points of a data set """
        chunk_sz = 1000
        chgpts = []
        # Iterate over chunks
        for i in range(0, len(times), chunk_sz):
                t = times[i:i+chunk_sz]
                N = len(t)
                T = t[-1]
                L0 = compute_likelihood(t)

                k = np.indices(t.shape)[0]
                Vk = t / T
                mu_k = -sum(1. / np.arange(k, N-1))
                mu_nk = -sum(1. / np.arange(N-k, N-1))
                v2_k = -sum(1. / np.arange(k, N-1)**2)
                v2_nk = -sum(1. / np.arange(N-k, N-1)**2)
                xi = pi**2/6 - sum(1. / np.arange(1, N-1)**2)
                sigma2_k = 4*k**2*v2_k + 4*(N-k)**2*v2_nk - 8*k*(N-k)*xi

                L0_E = -2*k*log(Vk) + 2*k*mu_k - 2*(N-k)*log(1-Vk) + 2*(N-k)*mu_nk
                Lbar = L0_E / sqrt(sigma2_k)
                W = 0.5*log(4*k*(N-k)/N**2)
                L = (L0 - mean(L0)) / + W

                k_max = argmax(L)
                Z = L[k_max]
                chgpts.append(k_max)


def test_discrete_data():
        times = []
        inten = []
        t = 0
        intensities = [ 1, 2, 3, 4, 3, 2, 1 ]
        for intensity in intensities:
                for n in range(1e3):
                        dt = np.random.exponential(1./intensity)
                        t += dt
                        inten.append(intensity)
                        times.append(t)
        return np.core.records.fromarrays([times, inten], names='t,I')

def test_gaussian_data(sigma=10, k=40, I0=10, I1=1000):
        times = []
        inten = []
        t = 0
        for n in range(2e4):
                I = I0 + I1 * np.exp(-(t-k)**2 / sigma**2)
                dt = np.random.exponential(1./I)
                t += dt
                inten.append(I)
                times.append(t)
        return np.core.records.fromarrays([times, inten], names='t,I')

if __name__ == '__main__':
        import sys
        #data = test_discrete_data()
        #data = test_gaussian_data()
        #times = data['t']

        times = np.fromfile(sys.argv[1], dtype=np.uint64)
        bins = np.core.records.fromfile(open(sys.argv[2]), formats='u8,u2', names='t,n')
        times = 1.*(times - times[0])

        L = compute_likelihood(times)

        import matplotlib.pyplot as pl
        ax = pl.subplot(211)
        ax.autoscale_view(scalex=True, tight=True)
        np.savetxt('hi', L)
        ax.plot(times, L)

        ax = pl.subplot(212)
        ax.plot(bins['t'], bins['n'], '+')
        ax.autoscale_view(scalex=True, tight=True)
        pl.show()

