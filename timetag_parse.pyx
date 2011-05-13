import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from cpython cimport bool

ctypedef unsigned char uint8_t
ctypedef unsigned long long uint64_t

DEF ENDIANESS="little"

IF ENDIANESS == "big":
        cdef uint64_t swap_record(uint64_t rec):
                return rec
ELIF ENDIANESS == "little":
        cdef uint64_t swap_record(uint64_t rec):
                cdef uint64_t ret
                cdef uint8_t* a = <uint8_t*> &rec
                cdef uint8_t* b = <uint8_t*> &ret
                b[0] = a[5]
                b[1] = a[4]
                b[2] = a[3]
                b[3] = a[2]
                b[4] = a[1]
                b[5] = a[0]
                return ret
ELSE:
        print 'Invalid ENDIANESS'

cdef packed struct StrobeEvent:
        np.uint64_t time
        np.uint8_t channels
strobe_event_dtype = np.dtype([('t', np.uint64), ('chs', np.uint8)])

def get_strobe_events(f, channel_mask, skip_wraps=0):
        cdef char* fname 
        if isinstance(f, str):
                fname = f
        else:
                fname = f.filename
        cdef FILE* fl = fopen(fname, "r")
        if fl == NULL:
                raise RuntimeError("Couldn't open file")

        if channel_mask > 0xf:
                raise RuntimeError("Invalid channel mask")
        cdef uint64_t mask = channel_mask << 36

        cdef size_t chunk_sz = 1024
        cdef unsigned int j = 0
        cdef uint64_t time_offset = 0

        cdef unsigned int wraps = 0
        cdef uint64_t rec
        cdef np.ndarray[StrobeEvent] chunk
        chunk = np.empty(chunk_sz, dtype=strobe_event_dtype)
        chunks = []

        while not feof(fl):
                res = fread(&rec, 6, 1, fl)
                if res != 1: break
                rec = swap_record(rec)

                # Handle timer wraparound
                wrapped = rec & (1ULL<<46) != 0
                wraps += wrapped
                if wraps < skip_wraps: continue
                if wrapped:
                        time_offset += (1ULL<<36)

                # Record event
                if rec & mask and not (rec & (1ULL << 45)):
                        t = rec & ((1ULL<<36)-1)
                        t += time_offset
                        chunk[j].time = t
                        chunk[j].channels = (rec >> 36) & 0xf
                        j += 1

                        # Start new chunk on filled
                        if j == chunk_sz:
                                chunks.append(chunk)
                                chunk = np.empty(chunk_sz, dtype=strobe_event_dtype)
                                j = 0
        
        chunks.append(chunk[:j])
        fclose(fl)
        return np.hstack(chunks)

cdef packed struct DeltaEvent:
        np.uint64_t time
        np.uint8_t channels
delta_event_dtype = np.dtype([('t', np.uint64), ('chs', np.uint8)])

def get_delta_events(f, skip_wraps=0):
        cdef char* fname 
        if isinstance(f, str):
                fname = f
        else:
                fname = f.filename
        cdef FILE* fl = fopen(fname, "r")
        if fl == NULL:
                raise RuntimeError("Couldn't open file")

        cdef size_t chunk_sz = 1024
        cdef unsigned int j = 0
        cdef uint64_t time_offset = 0

        cdef unsigned int wraps = 0
        cdef uint64_t rec
        cdef np.ndarray[DeltaEvent] chunk
        chunk = np.empty(chunk_sz, dtype=delta_event_dtype)
        chunks = [chunk]

        while not feof(fl):
                res = fread(&rec, 6, 1, fl)
                if res != 1: break
                rec = swap_record(rec)

                # Handle timer wraparound
                wrapped = rec & (1ULL<<46) != 0
                wraps += wrapped
                if wraps < skip_wraps: continue
                if wrapped:
                        time_offset += (1ULL<<36)

                # Record event
                if rec & (1ULL << 45):
                        t = rec & ((1ULL<<36)-1)
                        t += time_offset
                        chunk[j].time = t
                        chunk[j].channels = (rec>>36) & 0xf
                        j += 1

                        # Start new chunk on filled
                        if j == chunk_sz:
                                chunk = np.empty(chunk_sz, dtype=delta_event_dtype)
                                chunks.append(chunk)
                                j = 0
        
        fclose(fl)
        return np.hstack(chunks)

cdef packed struct DeltaSpan:
        uint64_t start_t, end_t
        uint8_t state
delta_span_dtype = np.dtype([
        ('start_t', np.uint64),
        ('end_t', np.uint64),
        ('state', np.uint8)])

def get_delta_spans(f, channel, skip_wraps=0):
        """
        Gets a list of delta channel spans. A span is defined as a period 
        during which the given delta channel has a static value.
        """
        cdef char* fname 
        if isinstance(f, str):
                fname = f
        else:
                fname = f.filename
        cdef FILE* fl = fopen(fname, "r")
        if fl == NULL:
                raise RuntimeError("Couldn't open file")

        cdef size_t chunk_sz = 1024
        cdef unsigned int j = 0
        cdef uint64_t time_offset = 0

        cdef uint64_t t
        cdef bool last_state = False
        cdef uint64_t last_t = 0

        cdef unsigned int wraps = 0
        cdef uint64_t rec
        cdef np.ndarray[DeltaSpan] chunk
        chunk = np.empty(chunk_sz, dtype=delta_span_dtype)
        chunks = [chunk]

        while not feof(fl):
                res = fread(&rec, 6, 1, fl)
                if res != 1: break
                rec = swap_record(rec)

                # Handle timer wraparound
                wrapped = rec & (1ULL<<46) != 0
                wraps += wrapped
                if wraps < skip_wraps: continue
                if wrapped:
                        time_offset += (1ULL<<36)

                # Record event
                state = ((rec>>(36+channel)) & 1) != 0
                if rec & (1ULL << 45) and state != last_state:
                        t = rec & ((1ULL<<36)-1)
                        t += time_offset
                        chunk[j].start_t = last_t
                        chunk[j].end_t = t
                        chunk[j].state = last_state
                        if last_t != 0: j += 1 # Throw out first span to get correct start time
                        last_t = t
                        last_state = state

                        # Start new chunk on filled
                        if j == chunk_sz:
                                chunk = np.empty(chunk_sz, dtype=delta_span_dtype)
                                chunks.append(chunk)
                                j = 0
        
        chunk[j].start_t = last_t
        chunk[j].end_t = t
        chunk[j].state = last_state
        j += 1

        fclose(fl)
        return np.hstack(chunks)

