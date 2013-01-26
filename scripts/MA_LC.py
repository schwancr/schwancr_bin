#!/home/schwancr/Installed/epd/bin/python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t',dest='Rdata_dir',help='Data directory to find XYZData for trajectories')
parser.add_option('-q',dest='Qdata_dir',help='Data directory to find QData for trajectories')
parser.add_option('--cr',dest='coef_rmsd',default=1,type=float, help='Coefficient for RMSD [ 1 ]')
parser.add_option('--cq',dest='coef_qnorm',default=1,type=float, help='Coefficient for Q-norm [ 1 ]')
parser.add_option('--gq',dest='Qgens_FN',default='./Gens_QData.npy',help='Generator\'s native contact array [ Gens_QData.npy ]')
parser.add_option('--gr',dest='Rgens_FN',default='./Gens_XYZData.lh5',help='Generator\'s XYZ trajectory [ Gens_XYZData.lh5 ]')
parser.add_option('-a',dest='ass_out',default='./Assignments.h5', help='Assignments output filename ( Can NOT exist now ). [ ./Assignments.h5 ]')
parser.add_option('-r',dest='dist_out',default='./Assignments.h5.Dist', help='Distance to generator output filename ( Can NOT exist now ). [ ./Assignments.h5.Dist ]')
parser.add_option('-N',dest='part',default=0,type=int,help='Partition if doing multiple assigning')
parser.add_option('-P',dest='procs',default=1,type=int,help='Number of procs used if doing multiple assigning')
parser.add_option('-i',dest='aind',help='Atom indices to use in calculating RMSDs')
options,args = parser.parse_args()

from numpy import *
from msmbuilder import Serializer, DistanceMetric, Trajectory
import os
import re
import multiprocessing as mp

RMSD = DistanceMetric.RMSD

QGens = load( options.Qgens_FN ).astype(uint8)
RGens = Trajectory.Trajectory.LoadFromLHDF( options.Rgens_FN )['XYZList']
if options.aind:
    Aind = loadtxt( options.aind )
    Rgens = Rgens[ :, Aind ]
else:
    Aind = None

print QGens.shape, RGens.shape

RGens = RMSD.PrepareData( RGens )
metrics = [ 'rmsd', 'qnorm' ]
coefficients = [ options.coef_rmsd, options.coef_qnorm ]


DistLC = DistanceMetric.LinearCombination( metrics, coefficients )
def AssignTraj( trajFNs ):
    print " Assigning trajectories: %s " % str( trajFNs )
    Rtrj = Trajectory.Trajectory.LoadFromLHDF( trajFNs[0] )['XYZList']
    if Aind: # If atom indices are given, then restrict to just those atoms.
        Rtrj = Rtrj[ :, Aind ]
    Rtrj = RMSD.PrepareData( Rtrj )
    Qtrj = load( trajFNs[1] ).astype( uint8 )
    N = len( Qtrj )
    Ass1d = []
    Dist1d = []
    for i in xrange( N ):
        Dists = DistLC.GetFastMultiDistance( [ Rtrj, Qtrj ], [ RGens, QGens ], i )
        Ass1d.append( Dists.argmin() )      
        Dist1d.append( Dists.min() )

    return Ass1d, Dist1d

def main():

    if options.procs == 1:
        Ass_outName = options.ass_out
        Dist_outName = options.dist_out
    else:
        Ass_outName = 'Assignments.%d.Ass' % options.part
        Dist_outName = 'Assignments.%d.RMSD' % options.part

    if os.path.exists( Ass_outName ):
        print "%s exists!" % Ass_outName
        exit()
    if os.path.exists( Dist_outName ):
        print "%s exists!" % Dist_outName

    QtrjList = [ ( int( file[3:-4] ), os.path.join( options.Qdata_dir,file )  ) for file in os.listdir( options.Qdata_dir ) if re.search( 'trj\d+\.npy', file )  ]
    QtrjList.sort()
    QtrjList = [ b for (a,b) in QtrjList ]

    RtrjList = [ ( int( file[3:-4] ), os.path.join( options.Rdata_dir,file )  ) for file in os.listdir( options.Rdata_dir ) if re.search( 'trj\d+\.lh5', file )  ]
    RtrjList.sort()
    RtrjList = [ b for (a,b) in RtrjList ]

    sol = []

    for i in range( len( RtrjList ) )[options.part::options.procs]: # Default is [0::1] so it will do all the data
        sol.append( AssignTraj( ( RtrjList[i], QtrjList[i] ) ) )

    Assignments = []
    Dists = []
    maxLength = 0

    for thing in sol:
        if len( thing[0] ) > maxLength:
            maxLength = len( thing[0] )
        Assignments.append( thing[0] )
        Dists.append( thing[1] )


    AssAry = -1 * ones( ( len( Assignments ), maxLength ) )
    DistAry = AssAry.copy()

    for i in range(len( Assignments ) ):
        AssAry[i][:len( Assignments[i] )] = Assignments[i] 
        DistAry[i][:len( Dists[i] )] = Dists[i]

    Serializer.SaveData( Ass_outName, AssAry.astype(int) )
    Serializer.SaveData( Dist_outName, DistAry )
    if options.procs != 1:
        Serializer.SaveData( 'Assignments.%d.WhichTrajs' % options.part, arange( len( RtrjList ) )[options.part::options.procs] )

if __name__ == '__main__':
    main()
