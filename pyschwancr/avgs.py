

import numpy as np

#@cython.boundscheck(False) # turn off bounds checking for entire function
def state_sums(flat_assigns, flat_data):

    assert len(flat_assigns) == len(flat_data)

    N = np.max(flat_assigns) + 1
    length = len(flat_assigns)
    running_sum = np.zeros(N, dtype=np.float)
    running_sum_sqr = np.zeros(N, dtype=np.float)

    for i in range(length):
        
        cur_state = flat_assigns[i]
        cur_data = flat_data[i]
        if cur_state < 0:
            continue

        running_sum[cur_state] += cur_data
        running_sum_sqr[cur_state] += cur_data * cur_data

    return running_sum, running_sum_sqr

    
