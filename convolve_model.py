from lmfit import Model
from scipy.signal import fftconvolve

class ConvolveModel(Model):
    def __init__(self, response, model, *args, **kwargs):
        def convolve(x, **kwargs):
            kwargs['x'] = x
            m = model.eval(**kwargs)
            a = fftconvolve(response, m, 'valid')
            a = a[:x.shape[0]]
            return a

        super(ConvolveModel, self).__init__(convolve, *args, **kwargs)
        # This is terrible
        self._func_allargs = model._func_allargs
        self._param_names = model.param_names

        self.response = response
        self.model = model
