#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-a',dest='ass_FN',default='./Assignments.h5',help='Assignments from msmbuilder [ ./Assignments.h5 ]')
parser.add_option('--gq',dest='Qgen_FN',default='./Gens_QData.npy',help='Gens with Qdata in them')
parser.add_option('--gr',dest='Rgen_FN',default='./Gens_XYZData.lh5',help='Gens with XYZ data in them')
parser.add_option('--sq',dest='Qnat_state',help='np.array with native contacts satisfied in the native state. If none then a vector of ones is generated.')
parser.add_option('--sr',dest='Rnat_state',help='PDB with the native state in it')
parser.add_option('--cr',dest='coef_rmsd',default=1,type=float,help='Coefficient for RMSD in the linear combination')
parser.add_option('--cq',dest='coef_qnorm',default=1,type=float,help='Coefficient for Q-Norm in the linear combination')
options, args = parser.parse_args()
 
from numpy import *
from msmbuilder import Serializer, DistanceMetric, Trajectory, Conformation
from pyschwancr import dataIO
import os, sys, re

# First load the trajectories
Ass = Serializer.LoadData( options.ass_FN ).astype(int)
Ass1d = Ass[ where( Ass >= 0 ) ].flatten()

if options.coef_qnorm != 0:
	QGens = dataIO.readData( options.Qgen_FN ).astype(uint8)
else:
	QGens = ones( ( Ass.max()+1, 1 ) ).astype(uint8)

RGens = Trajectory.Trajectory.LoadFromLHDF( options.Rgen_FN )

metrics = [ 'rmsd', 'qnorm' ]
coefficients = [ options.coef_rmsd, options.coef_qnorm ]


DistLC = DistanceMetric.LinearCombination( metrics, coefficients )

NatStateXYZ = Conformation.Conformation.LoadFromPDB( options.Rnat_state )

if options.Qnat_state:
	NatStateQData = dataIO.readData( options.Qnat_state )
else:
	NatStateQData = ones( QGens.shape[1] ).astype(uint8)

print QGens, NatStateQData

GenDists = DistLC.GetMultiDistance( [ RGens['XYZList'], QGens ], [ NatStateXYZ['XYZ'], NatStateQData ] )

GenDist2N = GenDists[ Ass1d ]

save('GenDist2N.npy', GenDist2N )
