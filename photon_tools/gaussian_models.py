from photon_tools.model_fit import Model, Parameter, register_model
from numpy import min, max, mean, exp, power, sqrt, log10, sum

@register_model('gaussian')
class GaussianModel(Model):
        """ Gaussian peak model. """
        params = [
                Parameter('offset',      'offset', def_value=0.0, def_scope='fitted'),
                Parameter('center',          'center', def_value=0.5, def_scope='fitted'),
                Parameter('width',          'width', def_value=0.5, def_scope='fitted'),
                Parameter('amplitude',      'amplitude', def_value=1.0, def_scope='fitted'),
        ]

        def __call__(self, p, x):
                offset = p['offset']
                center = p['center']
		width = p['width']
		amplitude = p['amplitude']

		return offset + amplitude * exp(-(x - center)**2 / (2*(width)**2))
		
@register_model('double_gaussian')
class DoubleGaussianModel(Model):
        """ Double Gaussian peak model. """
        params = [
                Parameter('offset',      'offset', def_value=0.0, def_scope='fitted'),
                Parameter('a_center',    'a_center', def_value=0.5, def_scope='fitted'),
                Parameter('a_width',     'a_width', def_value=0.5, def_scope='fitted'),
                Parameter('a_amplitude', 'a_amplitude', def_value=1.0, def_scope='fitted'),
		Parameter('b_center',    'b_center', def_value=0.5, def_scope='fitted'),
                Parameter('b_width',     'b_width', def_value=0.5, def_scope='fitted'),
                Parameter('b_amplitude', 'b_amplitude', def_value=1.0, def_scope='fitted'),
        ]

        def __call__(self, p, x):
                offset = p['offset']
                a_center = p['a_center']
		a_width = p['a_width']
		a_amplitude = p['a_amplitude']
                b_center = p['b_center']
		b_width = p['b_width']
		b_amplitude = p['b_amplitude']

		a_contribution = a_amplitude * exp(-(x - a_center)**2 / (2*(a_width)**2))
		b_contribution = b_amplitude * exp(-(x - b_center)**2 / (2*(b_width)**2))
	

		return offset + a_contribution + b_contribution 
