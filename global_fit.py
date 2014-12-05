from __future__ import division
from copy import deepcopy
import lmfit

class Curve(object):
    """
    This represents a curve/model pair within a GlobalFit.
    """
    def __init__(self, name, data, model, weights=None):
        self.data = data
        self.model = model
        self.name = name
        self.weights = weights

    def _residual(self, params, **kwargs):
        return self.model._residual(params, self.data, self.weights, **kwargs)

class GlobalFit(lmfit.Minimizer):
    """
    This represents a fit over a number of curves, each against its
    own model. These models can potentially share parameters (known as
    "tying" of the parameters).

    The parameters object accepted by global fits contain two types of
    parameters: per-curve parameters are parameters which correspond
    to a parameter of a curve, and global parameters which tie
    together several per-curve parameters. Per-curve parameters have
    the name of their curve prepended to their name to ensure
    uniqueness.
    """

    def __init__(self, curves={}, method='leastsq'):
        self.method = method
        self.curves = curves
        self.best_fit = None
        self.params = lmfit.Parameters()
        self.init_params = deepcopy(self.params)
        # A map from global parameter name to the sub-model that are tied to it
        self.global_params = {}
        lmfit.Minimizer.__init__(self, self._residual, self.params)

    def tie_all(self, param):
        """ Tie together all of the parameters of the given name. """
        tied = ['%s_%s' % (curve.name, param)
                for curve in self.curves.values()
                for name in curve.model.param_names
                if name == param]
        if len(tied) == 0:
            raise RuntimeError('No known curve has a parameter named %s' % param)
        self.global_params.setdefault(param, set()).update(tied)

    def tie(self, param, curves):
        """
        Tie together the named parameter for all of the listed curves.
        """
        self.global_params.setdefault(param, set()).update(tied)
        for curve in curves:
            if curve not in self.curves:
                raise RuntimeError('No curve named %s' % curve)
            tied = ['%s_%s' % (curve.name, param)
                    for name in curve.model.param_names
                    if name == param]
            if len(tied) == 0:
                raise RuntimeError('Curve %s does not have a parameter named %s' % (curve, param))
            self.global_params.add(tied[0])

    def add_curve(self, name, model, data, weights=None):
        """
        Add a curve to the fit.
        """
        if name in self.curves:
            raise RuntimeError('Curve named %s already exists' % name)
        curve = Curve(name, data, model, weights)
        self.curves[name] = curve
        
    def _all_tied(self):
        """
        Generate a dictionary of all tied submodel parameters,
        mapping each to its corresponding global parameters.
        """
        tied = {}
        for global_param, params in self.global_params.items():
            for p in params:
                tied[p] = global_param
        return tied

    def make_params(self):
        """
        All parameter tyings should be declared before this is called.

        The resulting Parameters object excludes tied parameters as these
        aren't free to be set.
        """
        params = lmfit.Parameters()
        tied = self._all_tied()
        # Add sub-model parameters
        for curve in self.curves.values():
            for name, param in curve.model.make_params().items():
                pname = '%s_%s' % (curve.name, name)
                if pname not in tied:
                    params[pname] = param

        # Add global parameters
        for name, tied in self.global_params.items():
            params[name] = lmfit.Parameter(name=name)

        return params

    def _real_params(self, params):
        """
        Generate a Parameters object for each curve
        """
        real = {}
        tied = self._all_tied()
        for curve in self.curves.values():
            p = curve.model.make_params()
            for name in p:
                pname ='%s_%s' % (curve.name, name)
                src = tied.get(pname, pname)
                p[name].value = params[src].value
            real[curve.name] = p
        return real

    def eval(self, params, **kwargs):
        """
        Return a dictionary containing the model of each curve under
        the given set of parameters.
        """
        fits = {}
        real = self._real_params(params)
        fits = {name: curve.model.eval(params=real[name], **kwargs)
                for name,curve in self.curves.items()}
        return fits

    def _residual(self, params, **kwargs):
        residuals = {}
        real = self._real_params(params)
        for curve in self.curves.values():
            r = curve.model._residual(real[curve.name],
                                      curve.data, curve.weights, **kwargs)
            residuals[curve.name] = r
        return np.hstack(residuals.values())
            
    def fit(self, params=None, *args, **kwargs):
        """
        Perform a fit.
        """
        if params is not None:
            self.params = params

        self.fcn_args = args
        self.userkws = kwargs
        self.init_params = deepcopy(self.params)
        self.init_fits = self.eval(params=self.init_params, **self.userkws)

        # run fit
        self.minimize() #method=self.method)
        self.best_fit = self.eval(params=self.params, **self.userkws)


test = True
if test:
    import numpy as np
    from lmfit.models import ExponentialModel
    xs = np.arange(1000)
    ys =  np.exp(-xs / 800)
    ys1 = 3*np.exp(-xs / 400) + ys
    ys2 = 10*np.exp(-xs / 100) + ys
    
    gf = GlobalFit()
    m = ExponentialModel(prefix='m_') # common model
    m1 = ExponentialModel(prefix='m1_') + m
    m2 = ExponentialModel(prefix='m2_') + m
    gf.add_curve('c1', m1, ys1)
    gf.add_curve('c2', m2, ys2)
    gf.tie_all('m_decay')
    params = gf.make_params()
    params['m_decay'].value = 100
    gf.fit(params, x=xs)
    print gf.params
    
    from matplotlib import pyplot as pl

    pl.plot(xs, ys1, label='curve 1')
    pl.plot(xs, ys2, label='curve 2')
    pl.plot(xs, gf.best_fit['c1'], label='fit 1')
    pl.plot(xs, gf.best_fit['c2'], label='fit 2')
    pl.legend()
    pl.show()
