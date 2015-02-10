from __future__ import division
import numpy as np
from numpy.linalg import norm
import scipy.optimize
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

    # Compute H assuming sigma=1
    Hno_scale = np.empty((Nmodels, Nmodels), dtype='f8')
    for n in range(Nmodels):
        for m in range(Nmodels):
            Hno_scale[m,n] = 2 / Npts * np.sum(models[m,:] * models[n,:])
    print np.linalg.eigvals(Hno_scale)

    # Pre-condition data
    Hinit = Hno_scale + np.diag(np.full(Nmodels,  1e-8))
    print np.linalg.cond(Hinit)
    b = np.dot(models, y)
    pChi = np.dot(np.linalg.inv(Hinit), b)
    yPred = np.dot(models.T, pChi)
    eta = yPred - y
    scale = 1e6 / (np.max(yPred) - np.min(yPred))
    eta = eta * scale
    yScaled = y * scale
    scaledVar = np.var(eta)

    # Introduce variance into H and g0
    H = np.empty((Nmodels, Nmodels), dtype='f8')
    g0 = np.empty(Nmodels, dtype='f8')
    for n in range(Nmodels):
        g0[n] = 2 / Npts * np.sum(models[n,:] * yScaled / scaledVar)
        for m in range(Nmodels):
            H[m,n] = 2 / Npts * np.sum(models[m,:] * models[n,:] / scaledVar)

    delta = lambda p: np.diag(np.log(p / expected) / p)
    q = lambda p: np.dot(g0, p) - np.dot(p, np.dot(H + nu * delta(p), p)) / 2

    for i in range(4):
        #print p
        delta = np.diag(np.log(p / expected) / p)
        q = lambda p: np.dot(g0, p) - np.dot(p, np.dot(H + nu * delta, p)) / 2

        if False:
            import matplotlib.pyplot as pl
            xs = np.linspace(1e-5, 1)
            for i in range(Nmodels):
                ys = []
                for x in xs:
                    p[i] = x
                    ys.append(q(p))
                pl.plot(xs, ys, label='%d'%i)

            pl.legend()
            pl.show()

        # Shrager's method
        #print H, delta
        #print 'cond', np.linalg.cond(H), np.linalg.cond(H + nu * delta)
        #pNew,_ = shrager(Q=(H + nu * delta), g=g0, C=-np.eye(Nmodels), d=np.zeros(Nmodels), x0=p, mu=0)
        #p = pNew

        #print 'rank', Nmodels, np.linalg.matrix_rank(np.hstack([H + nu * delta, -np.eye(Nmodels)]))
        if False:
            # CVXOPT quadratic programming
            import cvxopt.solvers
            kwargs = {'P': cvxopt.matrix(H + nu * delta),
                      'q': cvxopt.matrix(-g0),
                      'G': cvxopt.matrix(-np.eye(Nmodels)),
                      'h': cvxopt.matrix(np.zeros(Nmodels)),
                      'kktsolver': 'chol',
                      'initvals': p
                      #A=np.zeros((Nmodels,Nmodels)), b=np.zeros(Nmodels))['x']
            }
            res = cvxopt.solvers.qp(**kwargs)
            p = np.array(res['x'])[:,0]

        if False:
            # CVXOPT SOCP
            import cvxopt.solvers
            kwargs = {'P': cvxopt.matrix(H + nu * delta),
                      'q': cvxopt.matrix(-g0),
                      'G': cvxopt.matrix(-np.eye(Nmodels)),
                      'h': cvxopt.matrix(np.zeros(Nmodels)),
                      'initvals': p
                      #A=np.zeros((Nmodels,Nmodels)), b=np.zeros(Nmodels))['x']
            }
            res = cvxopt.solvers.qp(**kwargs)
            p = np.array(res['x'])[:,0]

        print p
        p = np.maximum(0, p)

        # g0 is gradient of chi^2
        sGrad = np.log(p / expected) + 1
        antipar = norm(sGrad / norm(sGrad) - g0 / norm(g0)) / 2
        #print p, q(p), antipar

        converged = antipar < delta_thresh
        if converged:
            return p

    return p

def simple_mem(y, models, sigma, p0=None, expected=None, nu=5e-6, delta_thresh=1e-4):
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

    # Compute hessian and gradient around 0
    H = np.empty((Nmodels, Nmodels), dtype='f8')
    g0 = np.empty(Nmodels, dtype='f8')
    for n in range(Nmodels):
        g0[n] = 2 / Npts * np.sum(models[n,:] * y / sigma**2)
        for m in range(Nmodels):
            H[m,n] = 2 / Npts * np.sum(models[m,:] * models[n,:] / sigma**2)

    def objective(p):
        delta = np.diag(nu * np.log(p / expected) / p)
        q = np.dot(g0, p) - np.dot(p, np.dot(H + delta, p)) / 2
        grad = g0 - np.dot(H + delta, p)
        return (-q, -grad)

    res = scipy.optimize.minimize(objective, p,
                                  method='TNC',
                                  bounds=[(1e-14, None) for p in range(Nmodels)],
                                  jac = True)
    return res.x
