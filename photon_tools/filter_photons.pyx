from libc.stdlib cimport *
from libc.stdio cimport *
from cpython cimport bool
from timetag_types cimport *
from timetag_types import *

def filter_by_spans(np.ndarray[StrobeEvent] strobes, np.ndarray[DeltaEvent] deltas):
        cdef size_t chunk_sz = 1024
        cdef unsigned int j = 0
        cdef uint64_t t_off = 0
        cdef np.ndarray[StrobeEvent] chunk
        chunk = np.empty(chunk_sz, dtype=strobe_event_dtype)
        chunks = []

        cdef unsigned int delta_idx = 0
        # HACK? Move to first span in which we are interested
        while deltas[delta_idx].state != True:
                delta_idx += 1
                if delta_idx == len(deltas):
                        raise RuntimeError('No deltas events with state==True')
        for i in range(strobes.shape[0]):
                while strobes[i].time >= deltas[delta_idx+1].start_t:
                        if delta_idx == len(deltas) - 2:
                                chunks.append(chunk[:j])
                                return np.hstack(chunks)

                        if deltas[delta_idx+1].start_t < deltas[delta_idx].start_t:
                                raise RuntimeError, "Data inconsistency: Decreasing delta timestamps (index %d: %d -> %d)" % \
                                        (delta_idx, deltas[delta_idx].start_t, deltas[delta_idx+1].start_t)

                        delta_idx += 1

                        # Compensate for discarded duration
                        if deltas[delta_idx].state == False:
                                # Assumes that states alternate
                                t_off += deltas[delta_idx+1].start_t - deltas[delta_idx].start_t


                if strobes[i].time >= deltas[delta_idx].start_t and deltas[delta_idx].state:
                        chunk[j].time = strobes[i].time - deltas[delta_idx].start_t + t_off
                        chunk[j].channels = strobes[i].channels
                        j += 1

                        # Start new chunk on filled
                        if j == chunk_sz:
                                chunks.append(chunk)
                                chunk = np.empty(chunk_sz, dtype=strobe_event_dtype)
                                j = 0

        chunks.append(chunk[:j])
        return np.hstack(chunks)

