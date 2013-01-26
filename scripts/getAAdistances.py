#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t',dest='traj_dir',default='./Trajectories', help='Directory to find lh5 trajectories [ ./Trajectories ]')
parser.add_option('-w',dest='write_dir',default='./Trajectories_AA',help='Directory to write output arrays to [ ./Trajectories_AA ]')
parser.add_option('-s',dest='scheme',default='CA',help='Scheme to use in calculating the distances. One of { "CA", "closest", "closest-heavy" } [ CA ]')
options, args = parser.parse_args()
 
import numpy as np
from Emsmbuilder import metrics
from pyschwancr import dataIO
from msmbuilder import Trajectory
import os, sys, re
 
# Check if folder exists:

if os.path.exists( options.write_dir ):
   print "Write directory (%s) already exists! Exiting..." % options.write_dir
   exit()

os.makedirs( options.write_dir )

trajList = dataIO.getTrajList( options.traj_dir )

getRightCap = { 'ca' : 'CA', 'closest' : 'closest', 'closest-heavy' : 'closest-heavy' }

CC = metrics.ContinuousContact( scheme = getRightCap[ options.scheme.lower() ] )

for trjFN in trajList:
   print "Working on %s" % trjFN 
   traj = Trajectory.Trajectory.LoadFromLHDF( trjFN )
   trajOut = CC.prepare_trajectory( traj )
   outFN =  trjFN.split('/')[-1][:-4] + '.npy'
   np.save( os.path.join( options.write_dir, outFN ), trajOut )
   del trajOut, traj

print "Done! Output saved to %s" % options.write_dir 
