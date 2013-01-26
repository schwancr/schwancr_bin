#!/usr/bin/env python
import numpy as np
from msmbuilder import Serializer, Project
import os, sys, re
from schwancrtools import ArgLib_E

options, prep_metric = ArgLib_E.parse(['projectfn','pcaobject'],new_arglist=[('pc_n','-N','--n-pc','Which eigenvector to project onto.', 0), ('outFN','--out','--output','output filename. Will use np.save',None) ], metric_parsers=True)

print "Using %s to prepare the trajectory. This should NOT be a pca metric" % str( prep_metric.__repr__() )


Proj = Project.LoadFromHDF( options.projectfn )
pca = Serializer.LoadFromHDF( options.pcaobject )
which_vec = int( options.pc_n )

v = pca['vecs'][ :, which_vec ]

print v.shape
data = []

for i in xrange( Proj['NumTrajs'] ):
    print "Working on trajectory %d" % i
    traj = Proj.LoadTraj( i )
    ptraj = prep_metric.prepare_trajectory( traj )
    ptraj = ptraj.conj() * v
    ptraj = ptraj.sum(axis=1)
    data.extend( ptraj )

data = np.array(data)

np.save( options.outFN, data )
