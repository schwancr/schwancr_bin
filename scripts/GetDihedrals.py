#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t',dest='traj_dir',default='./Trajectories',help='Directory containing the trajectories in lh5 format')
parser.add_option('-w',dest='write_dir',default='./Trajectories_Dihedrals',help='Output directory to write the results. [ ./Trajectories_Dihedrals ]' )
parser.add_option('-P',dest='procs',type=int,default=1,help='Number of processes to run [ 1 ]')

options, args = parser.parse_args()
 
import numpy as np
from pyschwancr import dataIO, msmTools
from msmbuilder import Trajectory
import os, sys, re
import multiprocessing as mp

def AnalyzeTraj( trajFN ):
	traj = Trajectory.Trajectory.LoadFromLHDF( trajFN )
	print "Working on trajectory %s" % trajFN,
	Psis, Phis, Omegas = msmTools.calcDihedrals( traj, procs=options.procs )
	trajName = trajFN.split('/')[-1][:-4]
	print "Saved output to: ",
	np.save( os.path.join( options.write_dir, '%s_Psis.npy' % trajName ), Psis )
	print os.path.join( options.write_dir, '%s_Psis.npy' % trajName ),
	np.save( os.path.join( options.write_dir, '%s_Phis.npy' % trajName ), Phis )
	print os.path.join( options.write_dir, '%s_Phis.npy' % trajName ),
	np.save( os.path.join( options.write_dir, '%s_Omegas.npy' % trajName ), Omegas )
	print os.path.join( options.write_dir, '%s_Omegas.npy' % trajName )
	del Psis, Phis, Omegas, traj
	return

if not os.path.isdir( options.write_dir ):
	os.mkdir( options.write_dir )

# First get the trajectory list
TrajList = [ os.path.join( options.traj_dir, fn ) for fn in os.listdir( options.traj_dir ) if re.search( '^trj\d+\.lh5', fn ) ]

for trajFN in TrajList:
	AnalyzeTraj( trajFN )

print "Done. Output saved to %s" % options.write_dir
