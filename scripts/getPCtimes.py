#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-c',dest='cov_matFN',help='Covariance matrix as calculated from getCovMat.py')
parser.add_argument('-p',dest='pcaFN',help='PCA dictionary with vecs and vals')
parser.add_argument('-o',dest='outFN',help='Output filename to save the timescales (in frames) for each mode. These are sorted the same as the vals are in the input pca object')
parser.add_argument('--dt',dest='dt',type=int,help='Time given to TrainPCA_cross.py. In frames')
options = parser.parse_args()
 
from msmbuilder import Serializer
import numpy as np
import os, sys, re

cov_mat = Serializer.LoadData( options.cov_matFN )

pca = Serializer.LoadFromHDF( options.pcaFN )

vals = pca['vals']
vals[np.where(vals<10**-8)] = 10**-8
vecs = pca['vecs']

var0s = np.array([ np.dot( np.dot( vecs[:,i:i+1].T, cov_mat ), vecs[:,i:i+1] )[0,0] for i in xrange( vecs.shape[1] ) ])
print var0s
timescales = - options.dt / np.log( vals / var0s )
print timescales
np.savetxt(options.outFN, timescales)
