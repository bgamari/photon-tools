from __future__ import division
import numpy as np
from .shrager import shrager

def mem(y, models, sigma, p0=None, expected=None, nu=5e-6, delta_thresh=1e-4):
    """
    Compute the maximum entropy mixture of models fitting observations.

    Follows Vinogradov and Wilson. *Applied Spectroscopy.* Volume 54, Number 6 (2000)

    :type y: array of shape ``(Npts,)``
    :param y: Observations
    :type models: array of shape ``(Nmodels, Npts)``
    :param models: Models
    :type sigma: array of shape ``(Npts,)``
    :param sigma: Standard error of points
    :type p0: array of shape ``(Nmodels,)``, optional
    :param p0: Initial guess at mixture weights
    :type expected: array of shape ``(Nmodels,)``, optional
    :param expected: Expected weights :math:`M`
    :type nu: ``float``
    :param nu: Regularization parameter
    :type delta_thresh: ``float``, optional
    :param delta_thresh: Maximum allowable anti-parallelism between
        gradients of :math:`\chi^2` and :math:`S` for convergence.

    """
    (Nmodels, Npts) = models.shape
    assert y.shape == (Npts,)
    assert sigma.shape == (Npts,)

    if expected is None:
        expected = 1
    else:
        assert expected.shape == (Nmodels,)

    # Initial value of p
    if p0 is None:
        p0 = np.ones(Nmodels)
    assert p0.shape == (Nmodels,)
    p = p0

    # Compute H and g^0
    H = np.empty((Nmodels, Nmodels), dtype=float)
    g0 = np.empty(Nmodels)
    for n in range(Nmodels):
        g0[n] = 2 / Npts * np.sum(models[n,:] * y / sigma**2)
        for m in range(Nmodels):
            H[m,n] = 2 / Npts * np.sum(models[m,:] * models[n,:] / sigma**2)
            
    while True:
        delta = np.diag(np.log(p / expected) / p)
        #q = np.dot(g0, p) - np.dot(p, np.dot(H + v * delta))
        p,_ = shrager(Q=(H + nu * delta), g=g0, C=-np.eye(Nmodels), d=np.zeros(Nmodels), x0=p)

        converged= False
        if converged:
            break

    return p
