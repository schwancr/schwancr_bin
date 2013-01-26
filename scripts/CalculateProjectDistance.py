#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p',dest='proj_FN',default='./ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option('-m',dest='metric',default='rmsd',help='Metric to use to calculate distance. Should be one of: rmsd, bool_cm, cont_cm, dihedral, pca')
parser.add_option('-i',dest='restInd',default='restInd',help='Restrictive data file for the particular metric. RMSD -> AtomIndices, CM -> Contact List, etc.')
parser.add_option('-s',dest='pdbFN',help='PDB filename to calculate all the distances to')
parser.add_option('-o',dest='outFN',help='Output filename to save data in .h5 format')
options, args = parser.parse_args()
 
import numpy as np
from Emsmbuilder import Project, Trajectory, metrics, Serializer
from pyschwancr import dataIO
import os, sys, re
 
# Set up the metric:
if os.path.exists( options.outFN ):
   print "Filename exists (%s)" % options.outFN

if options.metric.lower() == 'rmsd':
   Dist = metrics.RMSD( atomindices = np.loadtxt( options.restInd, int ) )
elif options.metric.lower() == 'bool_cm':
   Dist = metrics.BooleanContact( contacts = np.loadtxt( options.restInd, int ) )
elif options.metric.lower() == 'cont_cm':
   Dist = metrics.ContinuousContact( contacts = np.loadtxt( options.restInd, int ) )
elif options.metric.lower() == 'dihedral':
   Dist = metrics.Dihedral( angles = open( options.restInd ).read().strip() )
else:
   print "Need to enter one of rmsd, bool_cm, cont_cm, or dihedral as a metric"
   exit()

Proj = Project.Project.LoadFromHDF( options.proj_FN )

DistOut = np.ones( ( Proj['TrajLengths'].shape[0], Proj['TrajLengths'].max() ) ) * -1

trajList = dataIO.getTrajList( Proj['TrajFilePath'] )
pdb = Trajectory.Trajectory.LoadTrajectoryFile( options.pdbFN )
pdb = Dist.prepare_trajectory( pdb )
print pdb.shape
for i, trajFN in enumerate(trajList):
   print "Working on %s" % trajFN
   traj = Trajectory.Trajectory.LoadFromLHDF( trajFN )
   traj = Dist.prepare_trajectory( traj )
   print traj.shape
   DistOut[i, : Proj['TrajLengths'][i] ] = Dist.one_to_all( pdb, traj, 0 )

   del traj

Serializer.SaveData(options.outFN, DistOut )
