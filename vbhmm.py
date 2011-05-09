#/usr/bin/python

import numpy as np
from numpy import exp, log, sum, ones
from scipy.special import gammaln, psi

def kl_dirichlet(P, Q):
	alphaP = sum(P)
        alphaQ = sum(Q)
        return gammaln(alphaP) - gammaln(alphaQ) \
                - sum(gammaln(P)-gammaln(Q)) \
                + sum((P-Q) * (psi(P)-psi(P)))

def vbmm(emission_seqs, n_states, model, max_iter=100, tol=1e-4):
        total_length = sum(map(len, emission_seqs))
        # Hyperparameters
        alpha_a, alpha_b, alpha_pi = 1,1,1
        # pseudo-counts
        ua = ones(1,n_states) * (alpha_a/K)
        ub = ones(1,n_obs) * (alpha_b/n_obs)
        upi = ones(1,n_states) * (alpha_pi/n_states)

        # Pick an HMM from the prior to initialize the counts
        wa = array(shape=(n_states,n_states))
        wb = array(shape=(n_stataes,))
        from scipy.random.mtrand import dirichlet
        for i in range(n_states):
                wa[k,:] = dirichlet(ua,1) * total_length
                wb[k,:] = dirichlet(ub,1) * total_length
        wpi = dirichlet(upi,1) * len(emission_seqs)

        Fold = -Inf; ntol = tol*len(emission_seqs)
        for i in range(max_iter):
                # M step
                Wa = wa + np.tile(ua, (K, 1))
                Wb = wb + np.tile(ub, (K, 1))
                Wpi = wpi + upi

                astar = exp( psi(Wa) - np.tile(psi(sum(Wa)), (1,K)) )
                bstar = exp( psi(Wb) - np.tile(psi(sum(Wb)), (1,L)) )
                pistar = exp( psi(Wpi - psi(sum(Wpi))) )

                alpha, beta, scale = forward_backward(astar, bstar, pistar, data)
                lnZv = sum(log(scale),1)

                # Compute F
                Fa[i] = 0; Fb[i] = 0; Fpi[i] = 0;
                for kk in range(K):
                        Fa[i] -= kl_dirichlet(Wa[kk,:], ua)
                        Fb[i] -= kl_dirichlet(Wb[kk,:], ub)

                Fpi[i] = -kl_dirichlet(Wpi, upi)
                F[i] = Fa[i] + Fb[i] + Fpi[i] + lnZ[i]
        
                if F[i] < F[i-1] - 1e-6:
                        raise RuntimeError('Violation')
                elif (F[i]-F[2]) < (1+ntol)*(F[i-1]-F[2]) or not finite(F[i]):
                        model.Wa = Wa
                        model.Wb = Wb
                        model.Wpi = Wpi
                        return model

