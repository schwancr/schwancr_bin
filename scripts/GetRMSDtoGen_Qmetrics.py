#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option( '-a',dest='ass_FN', default='./Assignments.h5', help='Assignments.h5 from msmbuilder [ ./Assignments.h5 ]')
parser.add_option( '-p',dest='proj_FN', default='../ProjectInfo.h5', help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option( '-g',dest='gens_FN', help='Gens.lh5 from msmbuilder. This script will build this using QDists if you do not alreay have it.')
parser.add_option( '-q',dest='qdist_FN', default='./Assignments.h5.QDists', help='Assignments.h5.QDists from KCentersQDist.py [ ./Assignments.h5.QDists ]' )
parser.add_option( '-P',dest='procs', type=int, default=1, help='Number of processors to use [ 1 ]')
options, args = parser.parse_args()

from numpy import *
from msmbuilder import Serializer, Project, DistanceMetric
import multiprocessing as mp
import os
RMSD = DistanceMetric.RMSD

def AnalyzeState( stateInd ):
	global RMSDs
	whichTuple = where( Ass == stateInd )
	which = array( whichTuple ).T
	stateTraj = Proj.GetConformations( which )
	
	RMSDs[ whichTuple ] = RMSD.GetMultiDistance( stateTraj['XYZList'], Gens['XYZList'][ stateInd ] )

	print "Finished state %d" % stateInd 
	return

def GenerateGens():
	GenWhich = where( QDists == 0 ) 
	if len( GenWhich[0] ) != Ass.max() + 1:
		print "Cannot locate generator indices!"
		print "GenInd = %d. Number of States = %d" % ( len( GenWhich[0] ), Ass.max() + 1 )
		exit()
	print "Getting generator coordinates."
	Gens = Proj.GetConformations( array( GenWhich ).T )
	GenInd = Ass[ GenWhich ]
	print "Rearranging the generator trajectory."
	Gens['XYZList'] = Gens['XYZList'][ GenInd ] # Need to re-order the generators by the actual state indices
	Gens.SaveToLHDF( 'Gens.lh5' )
	print "Saved Generator trajectory to Gens.lh5"

	return Gens

def main():
	global Ass
	global Proj
	global RMSDs
	global QDists
	global Gens
	Ass = Serializer.LoadData( options.ass_FN )
	Proj = Project.Project.LoadFromHDF( options.proj_FN )
	QDists = Serializer.LoadData( options.qdist_FN )
	RMSDs = ones( Ass.shape ) * -1
	print "Loaded Data"
	if os.path.exists( options.ass_FN+'.RMSD' ):
		print "%s.RMSD Already Exists! Why are you trying to make it again?!" % options.ass_FN
		exit()

	try: Gens = Project.Trajectory.Trajectory.LoadFromLHDF( options.gens_FN )
	except: Gens = GenerateGens()
	
	pool = mp.Pool( options.procs )
	result = pool.map_async( AnalyzeState, range( Ass.max() + 1 ) )
	result.wait()

	Serializer.SaveData( options.ass_FN+'.RMSD', RMSDs )
	print "Saved Data to %s.RMSD" % options.ass_FN 


if __name__ == '__main__':
	main()
