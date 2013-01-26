#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-g',dest='gens_FN',default='./Gens.npy',help='Generator\'s native contact array')
parser.add_option('-d',dest='data_dir',help='Data directory to find assignments for trajectories')
parser.add_option('-p',dest='procs',type=int, default=1, help='Number of processors to use in assigning')
parser.add_option('-a',dest='ass_out',default='./Assignments.h5', help='Assignments output filename ( Can NOT exist now ). [ ./Assignments.h5 ]')
parser.add_option('-r',dest='dist_out',default='./Assignments.h5.QNorm', help='Distance to generator output filename ( Can NOT exist now ). [ ./Assignments.h5.QNorms ]')
options,args = parser.parse_args()

from numpy import *
from msmbuilder import Serializer, DistanceMetric
import os
import re
import multiprocessing as mp

Gens = load( options.gens_FN ).astype(uint8)

QNorm = DistanceMetric.QNorm
def AssignTraj( trajFN ):
    print " Assigning trajectory: %s " % trajFN
    trj = load( trajFN ).astype( uint8 )
    N = len( trj )
    Ass1d = []
    Qdist1d = []
    for i in xrange( N ):
        Dists = QNorm.GetMultiDistance( Gens, trj[i] )
        Ass1d.append( Dists.argmin() )      
        Qdist1d.append( Dists.min() )

    return Ass1d, Qdist1d

def main():

    Ass_outName = options.ass_out
    QDist_outName = options.dist_out
    if os.path.exists( Ass_outName ):
        print "%s exists!" % Ass_outName
        exit()
    if os.path.exists( QDist_outName ):
        print "%s exists!" % QDist_outName

    trjList = [ ( int( file[3:-4] ), os.path.join( options.data_dir,file )  ) for file in os.listdir( options.data_dir ) if re.search( 'trj\d+\.npy', file )  ]
    trjList.sort()
    trjList = [ b for (a,b) in trjList ]

    pool = mp.Pool( options.procs )
    result = pool.map_async( AssignTraj, trjList )
    result.wait()
    sol = result.get()

    Assignments = []
    QDists = []
    maxLength = 0

    for thing in sol:
        if len( thing[0] ) > maxLength:
            maxLength = len( thing[0] )
        Assignments.append( thing[0] )
        QDists.append( thing[1] )


    AssAry = -1 * ones( ( len( Assignments ), maxLength ) )
    QDistAry = AssAry.copy()

    for i in range(len( Assignments ) ):
        AssAry[i][:len( Assignments[i] )] = Assignments[i] 
        QDistAry[i][:len( QDists[i] )] = QDists[i]

    Serializer.SaveData( Ass_outName, AssAry.astype(int) )
    Serializer.SaveData( QDist_outName, QDistAry )


if __name__ == '__main__':
    main()
