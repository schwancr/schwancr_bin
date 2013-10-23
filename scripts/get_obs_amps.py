#!/usr/bin/env python

from msmbuilder import arglib
import numpy as np
from msmbuilder import io
from msmbuilder import msm_analysis
from scipy.io import mmread
import logging
logger = logging.getLogger("msmbuilder.arglib")

def run(tProb, observable, init_pops=None, num_vecs=10, output='evec_amps.h5'):

    if init_pops is None:
        init_pops = np.ones(tProb.shape[0]).astype(float) / float(tProb.shape[0])

    else:
        init_pops = init_pops.astype(float) 
        init_pops /= init_pops.sum()

    assert (observable.shape[0] == init_pops.shape[0])
    assert (observable.shape[0] == tProb.shape[0])
    
    try:
        f = io.loadh('eigs%d.h5' % num_vecs)
        vals = f['vals']
        vecsL = f['vecs']
    except:
        vals, vecsL = msm_analysis.get_eigenvectors(tProb, num_vecs + 1, right=False)
        io.saveh('eigs%d.h5' % num_vecs, vals=vals, vecs=vecsL)

    equil = vecsL[:,0] / vecsL[:,0].sum()

    dyn_vecsL = vecsL[:, 1:]
    # normalize the left and right eigenvectors

    dyn_vecsL /= np.sqrt(np.sum(dyn_vecsL * dyn_vecsL / np.reshape(equil, (-1, 1)), axis=0))

    dyn_vecsR = dyn_vecsL / np.reshape(equil, (-1, 1))

    amps = dyn_vecsL.T.dot(observable) * dyn_vecsR.T.dot(init_pops)

    io.saveh(output, evals=vals[1:], amplitudes=amps)
    logger.info("saved output to %s" % output)


if __name__ == '__main__':

    parser = arglib.ArgumentParser()
    parser.add_argument('tprob', help='transition probability matrix')
    parser.add_argument('observable', help='observable for each state in MSM')
    parser.add_argument('num_vecs', type=int, default=10, help='number of eigenvectors to look at')
    parser.add_argument('output', default='evec_amps.h5', help='output filename (use msmbuilder.io.loadh to read)')
    parser.add_argument('init_pops', default=None, help='initial populations')

    args = parser.parse_args()

    arglib.die_if_path_exists(args.output)

    ext = args.observable.split('.')[-1]
    if ext == 'h5':
        try:
            observable = io.loadh(args.observable, 'arr_0')
        except:
            observable = io.loadh(args.observable, 'Data')
    elif ext == 'npy':
        observable = np.load(args.observable)
    else:
        try:
            observable = np.loadtxt(args.observable)
        except:
            raise Exception("cannot read %s" % args.observable) 

    try:
        init_pops = np.loadtxt(args.init_pops).flatten().astype(float)
        init_pops /= init_pops.sum()
    except:
        init_pops = None

    if len(observable.shape) > 1:
        logger.warn("observable is two-dimensional. proceeding with the SECOND column")

        observable = observable[:, 1]
    
    tProb = mmread(args.tprob)

    if tProb.shape[0] != observable.shape[0]:
        raise Exception("observable does not match size of tProb (%d vs %d)" % (observable.shape[0], tProb.shape[0]))
    
    run(tProb, observable, init_pops=init_pops, num_vecs=args.num_vecs, output=args.output)
