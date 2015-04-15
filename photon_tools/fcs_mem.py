from __future__ import division
import numpy as np
from numpy.linalg import norm
import scipy.optimize

def fcs_mem(y, models, sigma, p0=None, expected=None, nu=5e-6, delta_thresh=1e-4):
    """
    Compute the maximum entropy mixture of models fitting observations.

    Very roughly follows Vinogradov and Wilson. *Applied
    Spectroscopy.* Volume 54, Number 6 (2000)

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
    :rtype: array of shape ``(Nmodels,)``
    :returns: Model weights
    """
    (Nmodels, Npts) = models.shape  # (N, M)
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
    p = p0.copy()

    # Compute hessian and gradient of \chi^2 at zero. Since we are working with
    # an additive mixture of components \chi^2 will contain at most second-order
    # terms in the weights. This implies that this expansion holds for
    # all points in the weight space.
    H = np.empty((Nmodels, Nmodels), dtype='f8')
    g0 = np.empty(Nmodels, dtype='f8')
    for n in range(Nmodels):
        g0[n] = 2 / Npts * np.sum(models[n,:] * y / sigma**2)
        for m in range(Nmodels):
            H[m,n] = 2 / Npts * np.sum(models[m,:] * models[n,:] / sigma**2)

    def objective_entropy(p):
        delta = np.diag(nu * np.log(p / expected) / p)
        q = np.dot(g0, p) - np.dot(p, np.dot(H + delta, p)) / 2
        grad = g0 - np.dot(H + delta, p)
        return (-q, -grad)

    def objective_l2(p):
        delta = np.diag(nu)
        q = np.dot(g0, p) - np.dot(p, np.dot(H + delta, p)) / 2
        grad = g0 - np.dot(H + delta, p)
        return (-q, -grad)

    res = scipy.optimize.minimize(objective_entropy, p,
                                  method='L-BFGS-B',
                                  bounds=[(1e-14, None) for p in range(Nmodels)],
                                  jac = True,
                                  tol = 1e-10)
    return res.x
