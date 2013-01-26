#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t',dest='traj_dir',default='./Trajectories',help='Directory to find trajectories [ ./Trajectories ]')
parser.add_option('-s',dest='pdbFN',help='PDB filename to align the trajectories to before doing PCA')
parser.add_option('-o',dest='outFN',default='pcaObj.npy',help='File to pickle the pca object to. This is an mdp.nodes.PCANode object. See mdp-toolkit for information about it. [ pcaObj.npy ]')
parser.add_option('-u',dest='stride',default=1,type=int,help='Stride to subsample at to train the PCA for. [ 1 ]')
parser.add_option('-N',dest='numComp',type=int,help='Number of components to find in PCA. NOTE: Only one of -N or -E should be set. Setting both will default to using -N"')
parser.add_option('-E',dest='explVar',type=float,help='Portion of explained variance to account for in picking how many dimensions to use from PCA. NOTE: Only one of -N or -E should be set. Setting both will default to using -N"')
options, args = parser.parse_args()

import numpy as np
from pyschwancr import ReduceDimensions
import os, sys, re
 
pca = ReduceDimensions.PCA(options.pdbFN, input_dir=options.traj_dir, expl_var=options.explVar, numComp=options.numComp, stride=options.stride)
print "Calculating the principal components"
pca.train()
pca.pca.stop_training()
print "PCA is finished. You kept %d components, which explains for a total of %f of the total variance" % (pca.pca.d.shape[0], pca.pca.explained_variance)
pca.savePCA(outFN=options.outFN)
