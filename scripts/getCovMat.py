#!/usr/bin/env python
 
from schwancrtools import ArgLib_E

options, prep_metric = ArgLib_E.parse(['projectfn'],new_arglist=[('outFN','-o','--out','outputfilename to save covariance matrix to. Will use Serializer.SaveData','./CovMat.h5')],metric_parsers=True)
 
from msmbuilder import Project, Serializer
import numpy as np
from schwancrtools import complexPCA
from pyschwancr import dataIO
import os, sys, re
 
Proj = Project.LoadFromHDF(options.projectfn)

cov_matrix = None
N = 0

for i in xrange(Proj['NumTrajs']):
    print "Working on trajectory %d" % i
    ptraj = prep_metric.prepare_trajectory( Proj.LoadTraj(i) ).astype(complex)

    n_rows, n_cols = ptraj.shape

    temp_cov_matrix = np.zeros( ( n_cols, n_cols ) ).astype(complex)

    if not ptraj.flags.contiguous:
        ptraj = ptraj.copy()
        if not ptraj.flags.contiguous:
            print "ptraj and ptraj.copy() are not contiguous arrays... Something wierd is going on."

    complexPCA.get_covariance( ptraj, temp_cov_matrix )
    temp_cov_matrix *= n_rows

    if cov_matrix == None:
        cov_matrix = temp_cov_matrix.copy()
    else:
        cov_matrix += temp_cov_matrix
   
    N += n_rows


cov_matrix /= N

Serializer.SaveData( options.outFN, cov_matrix )
print "Saved covariance matrix to %s" % options.outFN 
