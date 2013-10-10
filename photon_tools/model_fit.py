from collections import namedtuple
import numpy as np
from numpy import log10
import scipy
from scipy.optimize import leastsq

"""
model_fit - A flexible framework for model parameter fitting

This is a framework for doing non-linear least-squares fitting on several
datasets simultaneously.

"""

models = {}
def register_model(cls):
        if not issubclass(cls, Model):
                raise Exception("Registering model that doesn't inherit from model")
        models[cls.short_name] = cls
        return cls

class Parameter(object):
        def __init__(self, name, description, unit='', def_value=None, def_scope=None):
                self.name = name
                self.unit = unit
                self.description = description
                self.def_value = def_value
                self.def_scope = def_scope
                self.scope = None
                self.value = None

        def __str__(self):
                return '<Parameter %s (%s) = %e %s>' % (self.name, self.scope, self.value, self.unit)

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

        def __str__(self):
                return "\n".join(map(lambda x: str(x), self.values()))

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
        the __call__ function.
        """

        # Override this in subclasses
        params = []

        @classmethod
        def param(cls, name):
                for p in cls.params:
                        if p.name == name: return p
	
        def __call__(self, params, x):
                """ Compute the value of the fit function with the given
                    parameters on the given domain """
                pass

def _compute_error(p, curves, params, model):
        err = []
        params._unpack(p)
        for i,c in enumerate(curves):
                cparams = params._curve_params(i)
                G = model(cparams, c['lag'])
                err.extend((c['G'] - G) / np.sqrt(c['var']))

        err = np.array(err)
        return err

def fit(curves, model, params, epsfcn=0.0, verbose=False):
        """ Run the regression. Returns a new Model object with optimized
            parameters. One can then evaluate the optimized fit
            functions using this new object's G function. Curves is a
            list of record arrays. """
        from copy import deepcopy
        params = deepcopy(params)
        params.validate()
        p0 = params._pack()
        # TODO: Not sure why epsfcn needs to be changed
        res = scipy.optimize.leastsq(_compute_error, p0,
                                     args=(curves, params, model),
                                     full_output=True, epsfcn=epsfcn)
        p, cov_x, infodict, mesg, ier = res
        if verbose:
                print 'Completed after %d evaluations with message: %s' % (infodict['nfev'], mesg)
                print mesg
                print 'Error: ', sum(_compute_error(p, curves, params, model)**2)
                print ''
        if cov_x is None: raise RuntimeError('Fit failed to converge (flat axis)')
        params._unpack(p)
        return params, cov_x

def plot_model(fig, ax, params, model, curve_names, npts=1e3, with_uncertainty=False):
        """ plot_model(fig, ax, params, model, curve_names, npts=1e3, with_uncertainty=False)

            npts: Number of model points to plot
            uncertainty: Plot residuals with error bars showing uncertainty of data points
            curve_names: Curve labels for legend
        """
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        for i,curve in enumerate(params.curves):
                err = model(params._curve_params(i), curve['lag']) - curve['G']
                chi_sq = sum(err**2)
                dof = len(params.curves[i]) - len(params)

                name = curve_names[i].name
                start = log10(min(curve['lag']))
                stop = log10(max(curve['lag']))
                x = np.logspace(start, stop, npts)
                cparams = params._curve_params(i)
                m = model(cparams, x)
                ax.semilogx(x, m, label='%s (fit)' % name)

                ax2 = divider.append_axes('top', size=1.2, pad=0.1, sharex=ax)
                d = params.curves[i]
                resid = model(cparams, d['lag']) - d['G']
                if with_uncertainty:
                        ax2.errorbar(d['lag'], resid, marker='+', yerr=np.sqrt(d['var']))
                else:
                        ax2.scatter(d['lag'], resid, marker='+')
                ax2.set_xscale('log')
                # This was broken in matplotlib <1.2 due to Issue #1246
                ax2.axhline(y=0, color='k', alpha=0.5)
                for t in ax2.get_xticklabels():
                        t.set_visible(False)

                text =    'chi^2             = %1.1e' % chi_sq
                text += '\ndof               = %3d' % dof
                text += '\nchi^2/DOF         = %1.1e' % (chi_sq / dof)
                ptext = ["%3s  %-12s = %1.3e" % (p.scope[0:3], p.name, p.value)
                                for p in params.values() if isinstance(p.value, list)]
                text += '\n' + '\n'.join(sorted(ptext))
                ax2.text(1.05, 0.90, name,
                                weight='bold',
                                fontsize=8,
                                horizontalalignment='left',
                                verticalalignment='top',
                                transform=ax2.transAxes)
                ax2.text(1.07, 0.76, text,
                                fontsize=8,
                                family='monospace',
                                horizontalalignment='left',
                                verticalalignment='top',
                                transform=ax2.transAxes)

        ptext = ["%3s  %-12s = %1.3e" % (p.scope[0:3], p.name, p.value)
                        for p in params.values() if not isinstance(p.value, list)]
        text = '\n' + '\n'.join(sorted(ptext))
        ax.text(1.05, 0.95, 'Global Parameters',
                        weight='bold',
                        fontsize=8,
                        horizontalalignment='left',
                        verticalalignment='top',
                        transform=ax.transAxes)
        ax.text(1.07, 0.90, text,
                        fontsize=8,
                        family='monospace',
                        horizontalalignment='left',
                        verticalalignment='top',
                        transform=ax.transAxes)
