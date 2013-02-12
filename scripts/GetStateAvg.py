#!/usr/bin/env python


import numpy as np
from msmbuilder import io
import multiprocessing as mp
from msmbuilder import arglib
import time
import logging
logger = logging.getLogger('msmbuilder.scripts.GetStateAvg')
logger.setLevel(logging.DEBUG)

def get_state_avg(state_ind):
    """
    This function utilizes the GLOBABL variables ASSIGNMENTS_LIST and DATA_LIST
    defined in the boilerplate below
    """    
    #logger.info("Working on state %d", state_ind)
    state_data = np.concatenate([DATA_LIST[i][np.where(
                                 ASSIGNMENTS_LIST[i] == state_ind)]
                                 for i in xrange(len(DATA_LIST))])
    return (state_ind, state_data.mean(), state_data.std())

def run(output='stateAvg.dat', num_procs=1):

    # New function in this scope so that mp can access the assignment/data_lists.

    unique_states = np.unique(np.concatenate(
        [np.unique(assignments[np.where(assignments != -1)]) 
         for assignments in ASSIGNMENTS_LIST]))

    if num_procs == 1:  # Just do a for loop so errors are printed when debugging.
        out_data = np.array([get_state_avg(i) for i in unique_states])    
    else:
        pool = mp.Pool(num_procs)

        num_tasks = len(unique_states)
        result = pool.map_async(get_state_avg, unique_states)
     
        
        while not result._ready:
            logger.info("Working ... %d / %d tasks remaining", result._number_left * result._chunksize, num_tasks)
            time.sleep(10)

        out_data = np.array(result.get())

    np.savetxt(output, out_data)
    
    return

if __name__ == '__main__':

    parser = arglib.ArgumentParser()
    parser.add_argument('assignments', nargs='+', help='One or more assignments files')
    parser.add_argument('output')
    parser.add_argument('procs', type=int, default=1, 
                        help='Number of processors to use.')
    parser.add_argument('data', nargs='+', help='Raw data to average. There should be one data argument for each assignments argument.')

    args = parser.parse_args()

    if len(args.data) != len(args.assignments):
        raise Exception("Need to provide a data file for each assignments file.")

    ASSIGNMENTS_LIST = []
    for fn in args.assignments:
        try:
            assignments = io.loadh(fn)['arr_0']
        except:
            assignments = io.loadh(fn)['Data']

        ASSIGNMENTS_LIST.append(assignments)

    DATA_LIST = []
    for fn in args.data:
        try:
            data = io.loadh(fn)['arr_0']
        except:
            data = io.loadh(fn)['Data']

        DATA_LIST.append(data)

    for i in range(len(DATA_LIST)):
        if ASSIGNMENTS_LIST[i].shape != DATA_LIST[i].shape:
            raise Exception("Shape mismatch between %d'th assignments %s and data %s" % (i, str(ASSIGNMENTS_LIST[i].shape), str(DATA_LIST[i].shape)))

    run(output=args.output, num_procs=args.procs)

