#!/usr/bin/env python


import numpy as np
from argparse import ArgumentParser
from pyschwancr import averages
from msmbuilder import io
import os

def run(ass_fn_list, data_fn_list, out_fn):

    print ass_fn_list, data_fn_list, out_fn
    assert len(ass_fn_list) == len(data_fn_list)

    if os.path.exists(out_fn):
        raise Exception("path exists! (%s)" % out_fn)

    ass_list = []
    data_list = []
    for i in range(len(ass_fn_list)):
        ass_fn = ass_fn_list[i]
        data_fn = data_fn_list[i]
        try:
            ass = io.loadh(ass_fn)['arr_0'].astype(np.int)
        except:
            ass = io.loadh(ass_fn)['Data'].astype(np.int)

        try:
            data = io.loadh(data_fn)['arr_0'].astype(np.float)
        except:
            data = io.loadh(data_fn)['Data'].astype(np.float)

        if ass.shape != data.shape:
            print i, ass.shape.__repr__(), data.shape.__repr__()

        ass_1d = ass[np.where(ass != -1)]
        data_1d = data[np.where(ass != -1)]

        ass_list.append(ass_1d)
        data_list.append(data_1d)

    num_states = np.max([np.max(a) for a in ass_list]) + 1
    tot_sums = np.zeros(num_states)
    tot_sums_sqr = np.zeros(num_states)
    tot_pops = np.zeros(num_states)

    for i in range(len(ass_list)):
        cluster_pop = np.bincount(ass_list[i])

        sums, sums_sqr = averages.state_sums(ass_list[i], data_list[i])

        tot_sums[:len(sums)] += sums
        tot_sums_sqr[:len(sums_sqr)] += sums_sqr
        tot_pops[:len(cluster_pop)] += cluster_pop

    states = np.arange(num_states)    

    good_inds = np.where(tot_pops != 0)
    if good_inds[0].shape[0] != num_states:
        print "Some states have zero counts..."

    means = np.zeros(num_states)
    stds = np.zeros(num_states)

    means[good_inds] = tot_sums[good_inds] / tot_pops[good_inds].astype(np.float)
    stds[good_inds] = tot_sums_sqr[good_inds] / tot_pops[good_inds].astype(np.float) - np.power(means[good_inds], 2)
    stds = np.sqrt(stds)
    
    results = np.vstack((states, means, stds)).T

    np.savetxt(out_fn, results)

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('-a', dest='ass_fn', help='Assignments file', nargs='+')
    parser.add_argument('-o', dest='out_fn', help='Output filename')
    parser.add_argument('-d', dest='data_fn', help='Data filename for all conformations.', nargs='+')

    args = parser.parse_args()

    run(args.ass_fn, args.data_fn, args.out_fn)



