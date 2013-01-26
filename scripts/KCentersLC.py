#!/home/schwancr/Installed/epd/bin/python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-q',dest='Qdata_dir',help='Directory to look for numpy arrays of data (Q over time)')
parser.add_option('-t',dest='Tdata_dir',help='Directory to look for lh5\'s of XYZData over time')
parser.add_option('-u',type=int, default=1,dest='stride', help='Stride to sample by')
parser.add_option('-r',dest='cutoff',type=float, default=-1,help='Cutoff to quit clustering')
parser.add_option('-i',dest='aind_FN',default='../AtomIndices.dat', help='Atom index file to use in calculating RMSD')
parser.add_option('-k',dest='numGens',type=int,default=1000000,help='Maximum number of Generators')
parser.add_option('--cr',dest='coef_rmsd',type=float, default=1, help='Coefficient for rmsd [ 1 ]' )
parser.add_option('--cq',dest='coef_qnorm',type=float, default=1, help='Coefficeint for Q-norm [ 1 ]')

options, args = parser.parse_args()

from numpy import *
from msmbuilder import Clustering, Trajectory
import os
import re
from time import time
QList = []
XYZList = []
#from bitarray import bitarray
Aind = loadtxt(options.aind_FN,int)
QtrjList = [ ( int( file[3:-4] ), file ) for file in os.listdir( options.Qdata_dir ) if re.search( 'trj\d+\.npy', file )  ]
QtrjList.sort()
QtrjList = [ b for (a,b) in QtrjList ]

TtrjList = [ ( int( file[3:-4] ), file ) for file in os.listdir( options.Tdata_dir ) if re.search( 'trj\d+\.lh5', file )  ]
TtrjList.sort()
TtrjList = [ b for (a,b) in TtrjList ]

# The trj list is sorted now.
print "Loading QData Trajectories ..."
for trj in QtrjList:
#    print trj
    QList.extend( load( os.path.join( options.Qdata_dir, trj ) ).astype(uint8)[::options.stride] )

print "Loading XYZData Trajectories ..."
for trj in TtrjList:
    #XYZList.extend( Trajectory.Trajectory.LoadFromLHDF( os.path.join( options.Tdata_dir, trj ) )['XYZList'][::options.stride] )
    trjTemp = Trajectory.Trajectory.LoadFromLHDF( os.path.join( options.Tdata_dir, trj ) )
    trjTemp.RestrictAtomIndices( Aind )
    XYZList.extend( trjTemp['XYZList'] )

XYZList = array( XYZList )
QList = array( QList )

Metrics = [ 'rmsd', 'qnorm' ]
Coefficients = [ options.coef_rmsd, options.coef_qnorm ]
Datas = [ XYZList, QList ]

#QList = array([ bitarray(thing) for thing in QList ] )
print "K-Centers Clustering with %d conformations..." % QList.shape[0] 
a=time()
KCentersInd = Clustering.KCenters.ClusterLinearCombination( Datas, options.numGens, Metrics, Coefficients, Cutoff = options.cutoff )
b=time()
print b-a
save( 'Gens_QData.npy', QList[ KCentersInd ] )

GenTraj = Trajectory.Trajectory.LoadFromLHDF( os.path.join( options.Tdata_dir, TtrjList[0] ) )
GenTraj.RestrictAtomIndices( Aind )
GenTraj['XYZList'] = XYZList[ KCentersInd ] 
GenTraj.SaveToLHDF( 'Gens_XYZData.lh5' )
