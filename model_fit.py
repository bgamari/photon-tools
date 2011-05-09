from collections import namedtuple
import numpy as np
import scipy
from scipy.optimize import leastsq

models = {}
def register_model(*names):
        def reg(cls):
                if not issubclass(cls, Model):
                        raise Exception("Registering model that doesn't inherit from model")
                for n in names:
                        models[n] = cls
                return cls
        return reg

class Parameter(object):
        def __init__(self, name, description, def_value=None, def_scope=None):
                self.name = name
                self.description = description
                self.def_value = def_value
                self.def_scope = def_scope

class Parameters(dict):
        def __init__(self, model, curves):
                from copy import copy
                assert(len(curves) > 0)
                self.curves = curves
                
                for p in model.params:
                        a = copy(p)
                        a.scope = p.def_scope
                        a.value = p.def_value
                        self[p.name] = a

        def validate(self):
                """ Needs to be called before pack or unpack are used """
                self._fitted = []
                self._fixed = []
                for p in self.values():
                        if p.scope == 'fitted':
                                self._fitted.append(p.name)
                        elif p.scope == 'fixed':
                                self._fixed.append(p.name)
                        else:
                                raise RuntimeError('Invalid scope')

        def _unpack(self, packed):
                """ Unpack parameters from vector from leastsq """
                n = len(self.curves)
                i = 0
                for k in self._fitted:
                        if isinstance(self[k], list):
                                self[k].value = packed[i:i+n]
                                i += n
                        else:
                                self[k].value = packed[i]
                                i += 1

        def _pack(self):
                """ Pack fitted parameters into vector for leastsq """
                packed = []
                for k in self._fitted:
                        if isinstance(self[k].value, list):
                                packed.extend(self[k].value)
                        else:
                                packed.append(self[k].value)
                return np.array(packed)

        def _curve_params(self, curve):
                """ Generate full parameters dict from fixed and fitted parameters """
                params = {}
                for k in self:
                        v = self[k].value
                        params[k] = v[i] if isinstance(v, list) else v
                return params

class Model(object):
        """
        Represents a fitting model for a multi-curve non-linear regression
 
        A model has parameters, which can be either fixed or chosen for
        optimization. In the latter case, the parameter can be taken as
        determined by the fit, or fixed per-curve. In either case, the parameter
        can be set independently for each curve, or homogenous across all
        curves. Each parameter is set to be fixed/fitted and (in)homogenous at
        model creation.

        To implement a new fit function, one must inherit from the Model class,
        providing a class member 'params' list giving the parameters supported by
        the model and their default scope. The fit function itself is given by
        the compute_G function.
        """

        # Override this in subclasses
        params = []

        @classmethod
        def param(cls, name):
                for p in cls.params:
                        if p.name == name: return p
                
        def compute_G(self, params, x):
                """ Compute the value of the fit function with the given
                    parameters on the given domain """
                pass

def _compute_error(p, curves, params, model):
        err = []
        params._unpack(p)
        for i,c in enumerate(curves):
                cparams = params._curve_params(i)
                G = model.compute_G(cparams, c['lag'])
                #err.extend(c['G'] - G)
                err.extend((c['G'] - G) / c['var']**2)

        err = np.array(err)
        return err

def fit(curves, model, params):
        """ Run the regression. Returns a new Model object with optimized
            parameters. One can then evaluate the optimized fit
            functions using this new object's G function. """
        from copy import deepcopy
        params = deepcopy(params)
        params.validate()
        p0 = params._pack()
        # Not sure why epsfcn needs to be changed
        res = scipy.optimize.leastsq(_compute_error, p0, args=(curves, params, model), full_output=True, epsfcn=1e-7)
        p, cov_x, infodict, mesg, ier = res
        print cov_x, mesg
        #if cov_x is None: raise RuntimeError('Fit failed to converge (flat axis)')
        params._unpack(p)
        return params
