import sys
import numpy as np
cimport numpy as np

ctypedef unsigned long long uint64_t

def bin_photons(np.ndarray[np.uint64_t, ndim=1] times, int bin_width):
        cdef unsigned int chunk_sz = 10000
        cdef np.ndarray[np.uint16_t, ndim=1] chunk = np.empty(chunk_sz, dtype=np.uint16)
        chunks = []

        cdef Py_ssize_t i, j
        cdef uint64_t new_start
        cdef uint64_t bin_start = times[0]
        cdef short bin_count = 0
        cdef unsigned int bin = 0

        for i in range(times.shape[0]):
                if times[i] >= bin_start + bin_width:
                        new_start = (times[i] / bin_width) * bin_width
                        chunk[bin] = bin_count
                        bin += 1
                        if bin == chunk_sz:
                                chunks.append(chunk)
                                chunk = np.empty(chunk_sz, dtype=np.uint16)
                                bin = 0

                        # Account for zero bins
                        for j in range(bin_start+bin_width, new_start, bin_width):
                                chunk[bin] = 0
                                bin += 1
                                if bin == chunk_sz:
                                        chunks.append(chunk)
                                        chunk = np.empty(chunk_sz, dtype=np.uint16)
                                        bin = 0

                        bin_count = 0
                        bin_start = new_start

                bin_count += 1

        chunks.append(chunk[:bin])
        return np.hstack(chunks)

