import sys
import numpy as np
cimport numpy as np

ctypedef unsigned short uint16_t
ctypedef unsigned long long uint64_t

cdef packed struct Bin:
        uint64_t start_t
        uint16_t count
bin_dtype = np.dtype([
        ('start_t', np.uint64),
        ('count', np.uint16)])

def bin_photons(np.ndarray[np.uint64_t] times, int bin_width, bool include_zeros=True):
        cdef unsigned int chunk_sz = 10000
        cdef np.ndarray[Bin] chunk = np.empty(chunk_sz, dtype=bin_dtype)
        chunks = []

        cdef Py_ssize_t i, j
        cdef uint64_t new_start
        cdef uint64_t bin_start = times[0]
        cdef short bin_count = 0
        cdef unsigned int bin = 0

        for i in range(times.shape[0]):
                if times[i] >= bin_start + bin_width:
                        new_start = (times[i] / bin_width) * bin_width
                        chunk[bin].start_t = bin_start
                        chunk[bin].count = bin_count
                        bin += 1
                        if bin == chunk_sz:
                                chunks.append(chunk)
                                chunk = np.empty(chunk_sz, dtype=bin_dtype)
                                bin = 0

                        # Account for zero bins
                        if include_zeros:
                                for j in range(bin_start+bin_width, new_start, bin_width):
                                        chunk[bin].start_t = j
                                        chunk[bin].count = 0
                                        bin += 1
                                        if bin == chunk_sz:
                                                chunks.append(chunk)
                                                chunk = np.empty(chunk_sz, dtype=bin_dtype)
                                                bin = 0

                        bin_count = 0
                        bin_start = new_start

                bin_count += 1

        chunks.append(chunk[:bin])
        return np.hstack(chunks)

