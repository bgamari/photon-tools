# -*- coding: utf-8 -*- 

from __future__ import division
import numpy as np
from numpy import sqrt, min, max, sum, dot
from numpy.linalg import inv, norm

"""
Utilities to generate quadratic programs for testing
"""

def random_orthogonal(n, sigma=1):
    """
    Generate a random (in the Haar measure) orthogonal matrix by the naive QR method

    TODO: Rewrite by method of Stewart 1980.
    
    :type n: ``int``
    :param n: Dimension of matrix
    """
    x = np.random.normal(0, sigma, size=(n,n))
    q, r = np.linalg.qr(x, mode='complete')
    signs = np.diag(np.sign(r))
    return q * signs
    
    
class QuadProgram(object):
    pass

def random_quadratic_program(n_unknown, n_obs, cond_num, residual_ratio
                             n_eq,
                             n_active_start, n_active_sol, n_active_both, n_inactive):
    """
    Generate a random quadratic program by way of Lenard and Minkoff 1984.

    :type n_unknown: ``int``
    :param n_unknown: Number of unknown parameters
    :type residual_ratio: ``float``
    :param residual_ratio: Ratio of the final objective function to the initial objective function value (beta).
    """
    
    k0 = n_eq + n_active_start
    k1 = n_eq + n_active_sol + 1
    assert k0 <= n_unknown
    assert k1 <= n_unknown

    x0 = np.zeros(n_unknown)   # starting point
    xStar = np.ones(n_unknown) # optimal solution
    v = xStar - x0

    # Generate objective function
    p = random_orthogonal(n_obs)
    q = random_orthogonal(n_unknown)
    t = sqrt(cond_num)
    d = t**np.random.uniform(-1, +1, min(n_obs, n_unknown))
    D[0] = 1 / t
    D[-1] = t
    C = np.dot(p, D * q)

    def random_vec():
        a = np.random.normal(0, 1, n_unknown)
        return a / np.linalg.norm(a)
        
    eq_constrs = []
    for i in range(n_eq):
        a = random_vec()
        a = np.dot(np.eye(n_unknown) - np.outer(v,v) / dot(v,v), a)
        b = np.dot(a.T, x0)
        eq_constrs.append((a, '==', b))
        
    active_start_constrs = []
    for i in range(n_active_start):
        # Assume less than
        a = random_vec()
        b = np.dot(a, x0)
        s = '>=' if np.dot(a, xStar) >= b else '<='
        active_start_constrs.append((a, s, b))
        
    active_sol_constrs = []
    for i in range(n_active_sol):
        a = random_vec()
        b = np.dot(a, xStar)
        s = '>=' if np.dot(a, x0) >= b else '<='
        active_sol_constrs.append((a, s, b))

    active_both_constrs = []
    for i in range(n_active_both):
        a = random_vec()
        a = np.dot(np.eye(n_unknown) - np.outer(v,v) / dot(v,v), a)
        b = np.dot(a.T, x0)
        active_both_constrs.append((a, '<=', b))
        
    inactive_both_constrs = []
    for i in range(n_inactive_both):
        a = random_vec()
        u = np.random.uniform(0, 1)
        y = np.abs(np.random.normal())
        if np.random.uniform(0,1) <= 0.5:
            b = np.dot(a, x0) + y
            s = '<=' if np.dot(a, xStar) >= b else '>='
        else:
            b = np.dot(a, xStar) + y
            s = '<=' if np.dot(a, xStar) >= b else '>='
        inactive_both_constrs((a, s, b))

    # Generate residuals at optimium
    z = np.random.normal(size=n_unknown)
    eta = residual_ratio**2 / (1 - residual_ratio**2) / sum(z**2)
    tmp = dot(v, dot(C, z))
    betaSign = np.sign(dot(v, dot(C, z)))  # Ensure 〈v|h〉 > 0
    beta = eta * (tmp + betaSign * sqrt(tmp**2 + sum(np.dot(C,v)**2) / eta))
    rStar = beta * z  # residuals at solution
    d = dot(C, xStar) + rStar
    
    # Kuhn-Tucker constraint
    # All constraints satisfied by equality at solution
    J = eq_constrs + active_sol_constrs + active_both_constrs]
    lambdStar = [-1 if rho == '>=' else +1
                 for (a,rho,b) in J]
    g = dot(lambdStar, [a in for (a,_,_) in J])
    h = dot(C, rStar)
    assert dot(v, h) != 0

    # Scale factor
    if dot(v,g) > 0 and k1 > 1:
        alpha2 = k1  * dot(v, g) / (k1 - 1) / dot(v, h)
    else:
        alpha2 = 0.1 * norm(v) / norm(h)
        
    # TODO:
    #   Apply scale factor
    #   H
