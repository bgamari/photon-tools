# -*- coding: utf-8

from photon_tools.model_fit import Model, Parameter, register_model
from numpy import min, max, mean, exp, power, sqrt, log10, sum

@register_model('3d_diff')
class DiffusionModel(Model):
        """ Three-dimensional diffusion model. This includes the anamalous diffusion
            exponent alpha, which can be set to 1 for normal diffusion. """
        params = [
                Parameter('tau_d',      'Diffusion time', def_value=100, unit=u'μs', def_scope='fitted'),
                Parameter('a',          'Aspect ratio', def_value=3, def_scope='fitted'),
                Parameter('n',          'Concentration', def_value=0.5, def_scope='fitted'),
                Parameter('alpha',      'Anamalous diffusion exponent (1=normal diffusion)', def_value=1, def_scope='fixed'),
                Parameter('offset',     'Offset', def_value=0, def_scope='fixed'),
        ]

        def __call__(self, p, x):
                a = p['a']
                n = p['n']
                tau_taud = (x / (p['tau_d']*1e-6))**p['alpha']

                b = 1. / (1. + tau_taud)
                c = 1. / (1. + tau_taud * a**-2)
                return (1. / n) * b * sqrt(c) + p['offset']

@register_model('3d_diff_triplet')
class NormalDiffusionTripletModel(Model):
        """ Three-dimensional diffusion model with triplet correction """
        params = [
                Parameter('tau_d',      'Diffusion time constants', unit='μs', def_scope='fitted'),
                Parameter('tau_F',      'Triplet state relaxation time', unit='μs', def_scope='fitted'),
                Parameter('a',          'Aspect ratio', def_value=3, def_scope='fitted'),
                Parameter('n',          'Concentration', def_scope='fitted'),
                Parameter('F',          'Fraction of particles in triplet state'),
        ]

        def __call__(self, p, x):
                a = p['a']
                F = p['F']
                n = p['n']
                tau_taud = x / (p['tau_d'] * 1e-6)
                tau_tauF = x / (p['tau_F'] * 1e-6)

                b = 1. / (1. + tau_taud)
                c = 1. / (1. + tau_taud * a**-2)
                d = (1. - F + F * exp(-tau_tauF)) / (1. - F)
                return d * b * sqrt(c) / n

