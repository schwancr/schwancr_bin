

import numpy as np
cimport numpy as np
cimport cython

DTYPE_FLOAT = np.float
ctypedef np.int_t DTYPE_INT_t
ctypedef np.float_t DTYPE_FLOAT_t

@cython.boundscheck(False) # turn off bounds checking for entire function
def state_sums(np.ndarray[DTYPE_INT_t, ndim=1] flat_assigns, 
               np.ndarray[DTYPE_FLOAT_t, ndim=1] flat_data):

    assert len(flat_assigns) == len(flat_data)

    cdef int N = np.max(flat_assigns)
    cdef int length = len(flat_assigns)
    cdef unsigned int i
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] running_sum = np.zeros(N, dtype=DTYPE_FLOAT)
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] running_sum_sqr = np.zeros(N, dtype=DTYPE_FLOAT)
    cdef double cur_data
    cdef int cur_state

    for i in range(length):
        
        cur_state = flat_assigns[i]
        cur_data = flat_data[i]
        if cur_state < 0:
            continue

        running_sum[cur_state] += cur_data
        running_sum_sqr[cur_state] += cur_data

    return running_sum, running_sum_sqr

    
