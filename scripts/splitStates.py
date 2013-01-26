#!/home/schwancr/Installed/epd/bin/python -u
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option('-a',dest='ass_FN',default='./Assignments.h5',help='Assignments from msmbuilder [ ./Assignments.h5 ]')
parser.add_option('-g',dest='gens_FN',default='./Gens.lh5',help='Generators from msmbuilder [ ./Gens.lh5 ]')
parser.add_option('-d',dest='data_FN',help='Data to split states with')
parser.add_option('-c',dest='cutoff',default=0.5,type=float,help='Cutoff to split states on [ 0.5 ]' )
parser.add_option('-w',dest='write_dir',default='./SplitStates',help='Directory to write output ( new assignments, new generators, map file ) [ SplitStates ]' )
parser.add_option('-P',dest='procs',default=1,type=int,help='Number of processes to run [ 1 ]')
parser.add_option('-s',dest='statesFN',help='States to split ( if you already know them )' )
options, args = parser.parse_args()
 
import numpy as  np
import multiprocessing as mp
from msmbuilder import Serializer, Trajectory, Project, DistanceMetric
from pyschwancr import dataIO, msmTools
import os, sys, re

def AnalyzeState( state ):
   print "Working on state %d" % state
   stateInd = np.where( ass == state )
   if not belowCut[ stateInd ].sum() in [ 0, stateInd[0].shape[0] ]: # split the state
      return True
   else:
      return False

   return False


# Read in the data
proj = Project.Project.LoadFromHDF( options.proj_FN )
data = dataIO.readData( options.data_FN )
ass = Serializer.LoadData( options.ass_FN )
rmsd = Serializer.LoadData( options.ass_FN + '.RMSD' ) # This could end poorly.... but I can add another parameter if need be...
gens = Trajectory.Trajectory.LoadFromLHDF( options.gens_FN )
print "Loaded the data"
if os.path.isdir( options.write_dir ):
   print "Directory exists, will write data to %s... Careful, this could cause options, since overwriting will CRASH this script" % options.write_dir
else:
   os.mkdir( options.write_dir )
   print "Made output directory (%s)" % options.write_dir

data2d = msmTools.reshapeRawData( data, proj )

maxState = ass.max()

belowCut = ( data2d <= options.cutoff )
aboveCut = ( data2d > options.cutoff )

if not options.statesFN:
   pool = mp.Pool( options.procs )
   result = pool.map_async( AnalyzeState, range( maxState + 1 ) )
   sol = result.get()

   toSplit = np.where( np.array( sol ) )[0]
   np.savetxt( os.path.join( options.write_dir, 'SplitTheseStates.dat' ),toSplit,"%d" )
else:
   toSplit = dataIO.readData( options.statesFN ).astype(int)
   if toSplit.shape == ():
      toSplit = np.array( [ toSplit ] )
#print toSplit
print "Splitting %d states" % len( toSplit )
# Splitting states, strategy:
#  1) Find out where the generator lies. The new state will be defined by the opposite
#  2) Change the assignments for the split conformations
#  3) Get the new generator based on whatever was closest to the previous generator
#  4) Add the conformation to the gens trajectory
#  5) Get the RMSD to the generators from each of the conformations
#  6) Update the RMSD file with the results from 5

newStateNum = maxState + 1
map = []

for state in toSplit: 
   print "Working on state %d" % state
   map.append( ( state, newStateNum ) ) # map contains an array mapping which states were split/added
   stateInd = np.where( ass == state )
   genInd = np.where( ( abs( rmsd ) <= 0.01 ) * ( ass == state ) )
   if belowCut[ genInd ]: # If the generator is below the cutoff, then grab the conformations above it:
      newStateInd = np.where( ( aboveCut ) * ( ass == state ) ) # where the confs are above the cutoff and the assignment is this state
   else:
      newStateInd = np.where( ( belowCut ) * ( ass == state ) ) # Else where the confs are below the cutoff and the assignment is this state
   # Now rewrite the assignments
   ass[ newStateInd ] = newStateNum
   # Get the new generator from  the minimum rmsd to the previous generator
   minRMSDtoPrevGen = rmsd[ newStateInd ].min()
   newGenInd = np.where( ( rmsd <= ( minRMSDtoPrevGen + 0.01 ) ) * ( rmsd >= ( minRMSDtoPrevGen - 0.01 ) ) * ( ass == newStateNum ) ) # Where is the conformation with the lowest rmsd to the previous generator, and is in the newest state
   try: newGen = proj.GetConformations( np.array( [ np.array( newGenInd ).T[0] ] ) ) # Get the conformation's coordinates. In case multiple conformations are close, just pick the first one.
   except: 
      print "No generator found for state %d, continuing ..." % state
      map.pop( -1 ) # Skip this, so undo what we already did
      ass[ newStateInd ] = state
      continue
   gens['XYZList'] = np.concatenate( ( gens['XYZList'], newGen['XYZList'] ) ) # Add the new generator to the gens trajectory
   
   tempTraj = proj.GetConformations( np.array( newStateInd ).T ) # Get the conformations that we are splitting off
   newRMSDtoGen = DistanceMetric.RMSD.GetMultiDistance( tempTraj['XYZList'], newGen['XYZList'][0] ) # Get RMSD to the new generator for the split state
   rmsd[ newStateInd ] = newRMSDtoGen # save the new RMSDs to the rmsd to generator object
   
   newStateNum += 1 # next state! 


# Save the files.
print "Saving data to %s" % options.write_dir
gens.SaveToLHDF( os.path.join( options.write_dir, 'Gens.lh5.Split' ) )
Serializer.SaveData( os.path.join( options.write_dir, 'Assignments.h5.Split' ), ass )
Serializer.SaveData( os.path.join( options.write_dir, 'Assignments.h5.Split.RMSD' ), rmsd )
np.savetxt( os.path.join( options.write_dir, 'SplitMap.dat' ), np.array( map ), "%d" )
