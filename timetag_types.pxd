cimport numpy as np

ctypedef unsigned char uint8_t
ctypedef unsigned long long uint64_t

cdef packed struct DeltaEvent:
        uint64_t start_t
        uint8_t state

cdef packed struct StrobeEvent:
        np.uint64_t time
        np.uint8_t channels

