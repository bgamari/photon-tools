cimport libc.stdint
cimport numpy as np
import numpy as np

cdef extern from "picoharp_helper.h":
    timestamps picoharp_read_file(char* filename, unsigned int channel)

    struct timestamps:
            double jiffy
            np.uint64_t* timestamps
            unsigned int nrec

def read_timestamps(filename, channel):
    ts = picoharp_read_file(filename, channel)
    n = ts.nrec
    cdef data = np.asarray(<np.uint64_t[:n]> ts.timestamps)
    return (data, ts.jiffy)
