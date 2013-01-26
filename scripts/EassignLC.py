#!/home/schwancr/Installed/epd/bin/python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t',dest='Rdata_dir',help='Data directory to find XYZData for trajectories')
parser.add_option('--cr',dest='coef_rmsd',default=1,type=float, help='Coefficient for RMSD [ 1 ]')
parser.add_option('--cq',dest='coef_qnorm',default=1,type=float, help='Coefficient for Q-norm [ 1 ]')
parser.add_option('-g',dest='gens_FN',default='./gens.lh5',help='Generator\'s XYZ trajectory [ Gens_XYZData.lh5 ]')
parser.add_option('-a',dest='ass_out',default='./Assignments.h5', help='Assignments output filename ( Can NOT exist now ). [ ./Assignments.h5 ]')
parser.add_option('-r',dest='dist_out',default='./Assignments.h5.Dist', help='Distance to generator output filename ( Can NOT exist now ). [ ./Assignments.h5.Dist ]')
parser.add_option('--ai',dest='aind_FN',help='Atom indices to apply to raw trajectories.')
parser.add_option('--ri',dest='rind_FN',help='Residue indices to apply for Q-Norm calculation.')
options,args = parser.parse_args()

from numpy import *
from Emsmbuilder import Serializer, metrics, Trajectory
import os
import re
import multiprocessing as mp

Aind = loadtxt( options.aind_FN, int )
Rind = loadtxt( options.rind_FN, int )
Rind -= 1
Coefficients = [ options.coef_rmsd, options.coef_qnorm ]
RMSD = metrics.RMSD( atomindices=Aind )
QNorm = metrics.BooleanContact( contacts=Rind, cutoff=0.6 )
DistLC = metrics.Hybrid( [ RMSD, QNorm ], Coefficients  )

RGens = Trajectory.Trajectory.LoadFromLHDF( options.gens_FN )
RGens = DistLC.prepare_trajectory( RGens )

def AssignTraj( trajFNs ):
    print " Assigning trajectories: %s " % str( trajFNs )
    Rtrj = Trajectory.Trajectory.LoadFromLHDF( trajFNs )
    N = Rtrj['XYZList'].shape[0]

    Rtrj = DistLC.prepare_trajectory( Rtrj )

    Ass1d = []
    Dist1d = []
    for i in xrange( N ):
        Dists = DistLC.one_to_all( Rtrj, RGens, i )
        Ass1d.append( Dists.argmin() )      
        Dist1d.append( Dists.min() )

    return Ass1d, Dist1d

def main():

    Ass_outName = options.ass_out
    Dist_outName = options.dist_out
    if os.path.exists( Ass_outName ):
        print "%s exists!" % Ass_outName
        exit()
    if os.path.exists( Dist_outName ):
        print "%s exists!" % Dist_outName

    RtrjList = [ ( int( file[3:-4] ), os.path.join( options.Rdata_dir,file )  ) for file in os.listdir( options.Rdata_dir ) if re.search( 'trj\d+\.lh5', file )  ]
    RtrjList.sort()
    RtrjList = [ b for (a,b) in RtrjList ]

    sol = []

    for i in range( len( RtrjList ) ):
        sol.append( AssignTraj( RtrjList[i] ) )

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


if __name__ == '__main__':
    main()
