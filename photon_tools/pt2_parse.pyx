import numpy as np
cimport numpy as np
from cpython cimport PyObject, Py_INCREF
from libc.stdlib cimport free
from libc.stdint cimport uint64_t

np.import_array()

cdef extern from "pt2.h":
    uint64_t *get_pt2_timestamps(char*, unsigned int, unsigned int*)

cdef class ArrayWrapper:
    cdef void* data_ptr
    cdef int size

    cdef set_data(self, int size, void* data_ptr):
        """ Set the data of the array

	    This cannot be done in the constructor as it must recieve C-level
	    arguments.

	    Parameters:
	    -----------
	    size: int
	    Length of the array.
	    data_ptr: void*
	    Pointer to the data
	"""

        self.data_ptr = data_ptr
        self.size = size

    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
	    tries to get an array from the object."""
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape,
                                               np.NPY_UINT64, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
	    references to the object are gone. """
        free(<void*>self.data_ptr)

def get_photon_events(char* filename, unsigned int channel):
    cdef uint64_t *array
    cdef np.ndarray ndarray
    cdef unsigned int n
    # Call the C function
    array = get_pt2_timestamps(filename, channel, &n)
    size = 8 * n;

    array_wrapper = ArrayWrapper()
    array_wrapper.set_data(size, <void*> array)
    ndarray = np.array(array_wrapper, copy=False)
    # Assign our object to the 'base' of the ndarray object
    ndarray.base = <PyObject*> array_wrapper
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(array_wrapper)

    return ndarray

