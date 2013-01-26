#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-q',dest='Q_dir',default='../Trajectories_NC',help='Directory with np.arrays of Native contacts for each trajectory')
parser.add_option('-t',dest='XYZ_dir',default='../Trajectories',help='Directory with Trajectory objects saved as .lh5')
parser.add_option('--sq',dest='Qnat_state',help='np.array with native contacts satisfied in the native state. If none then a vector of ones is generated.')
parser.add_option('--sr',dest='Rnat_state',help='PDB with the native state in it')
parser.add_option('--cr',dest='coef_rmsd',default=1,type=float,help='Coefficient for RMSD in the linear combination')
parser.add_option('--cq',dest='coef_qnorm',default=1,type=float,help='Coefficient for Q-Norm in the linear combination')
parser.add_option('-o',dest='output',default='RawDist.npy',help='Output filename for data (np.save format) [ RawDist.npy ]')
options, args = parser.parse_args()
 
import numpy as np
from msmbuilder import Serializer, DistanceMetric, Trajectory, Conformation
from pyschwancr import dataIO
import os, sys, re

def AnalyzeTraj( TrajInd ):
	print "Working on Trajectories %s and %s" % (QTrajs[ TrajInd ], RTrajs[ TrajInd ] )
	QTraj = dataIO.readData( QTrajs[ TrajInd ] ).astype(np.uint8)
	RTraj = Trajectory.Trajectory.LoadFromLHDF( RTrajs[ TrajInd ] )
	
	TrajDist = DistLC.GetMultiDistance( [ RTraj['XYZList'], QTraj ], [ NatStateXYZ['XYZ'], NatStateQData ] )

	return TrajDist

# First create a list of the trajectories

RTrajs = dataIO.getTrajList( options.XYZ_dir )
QTrajs = dataIO.getTrajList( options.Q_dir, RegEx = r'^trj\d+\.npy' )

if len( RTrajs ) != len( QTrajs ):
	print "Need the same number of trajectories in XYZ_Dir (%s) and Q_Dir (%s)" % ( options.XYZ_dir, options.Q_dir )
	exit()

metrics = [ 'rmsd', 'qnorm' ]
coefficients = [ options.coef_rmsd, options.coef_qnorm ]
DistLC = DistanceMetric.LinearCombination( metrics, coefficients )

NatStateXYZ = Conformation.Conformation.LoadFromPDB( options.Rnat_state )
if options.Qnat_state:
	NatStateQData = dataIO.readData( options.Qnat_state )
else:
	# Need to get dimension of the QData to generate the native state:
	Qn = dataIO.readData( QTrajs[0] ).shape[1]
	NatStateQData = np.ones( Qn ).astype(np.uint8)

TotalDists = []

for i in range( len( RTrajs ) ):
	TotalDists.extend( AnalyzeTraj( i ) )

TotalDists = np.array( TotalDists )	

np.save( options.output, TotalDists )
