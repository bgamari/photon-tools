import numpy as np
cimport numpy as np

ctypedef unsigned long long uint64_t
def bin_photons(np.ndarray[np.uint64_t, ndim=1] times, int bin_width, bool no_zeros):
        cdef nbins = (times[-1] - times[0]) / bin_width + 1
        cdef np.ndarray[np.uint16_t, ndim=1] bins = np.empty(nbins, dtype=np.uint16)
        cdef Py_ssize_t i, j
        cdef uint64_t new_start
        cdef uint64_t bin_start = times[0]
        cdef short bin_count = 0
        cdef unsigned int bin = 0

        for i in range(times.shape[0]):
                if times[i] > bin_start + bin_width:
                        new_start = (times[i] / bin_width) * bin_width
                        if not no_zeros:
                                for j in range(bin_start+bin_width, new_start, bin_width):
                                        bins[bin] = 0
                                        bin += 1
                        
                        bins[bin] = bin_count
                        bin += 1
                        bin_count = 0
                        bin_start = new_start
                bin_count += 1

        return bins

