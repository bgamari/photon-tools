from lmfit import Model
from scipy.signal import fftconvolve

class ConvolveModel(Model):
    def __init__(self, response, model, *args, **kwargs):
        def convolve(x, **kwargs):
            #model_args = {k: v for k,v in kwargs.items() if k in model.param_names}
            model_args = model.make_funcargs(kwargs=kwargs)
            m = model.func(x, **model_args)
            a = fftconvolve(response, m, 'valid')
            a = a[:x.shape[0]]
            return a

        super(ConvolveModel, self).__init__(convolve, *args, **kwargs)
        # This is terrible
        self._func_allargs = model._func_allargs
        self._param_names = model.param_names

        self.response = response
        self.model = model
