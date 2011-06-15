import numpy as np

strobe_event_dtype = np.dtype([('t', np.uint64), ('chs', np.uint8)])
delta_event_dtype = np.dtype([
        ('start_t', np.uint64),
        ('state', np.uint8)])

bin_dtype = np.dtype([
        ('start_t', np.uint64),
        ('count', np.uint16)])
