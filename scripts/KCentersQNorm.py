#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='data_dir',help='Directory to look for numpy arrays of data (Q over time)')
parser.add_option('-u',type=int, default=1,dest='stride', help='Stride to sample by')
parser.add_option('-r',dest='cutoff',type=float, default=-1,help='Cutoff (Number of dissimilar native contacts)')
parser.add_option('-k',dest='numGens',type=int,default=1000000,help='Maximum number of generators to create')
options, args = parser.parse_args()

from numpy import *
from msmbuilder import Clustering
import os
import re
from time import time
QList = []
#from bitarray import bitarray

trjList = [ ( int( file[3:-4] ), file ) for file in os.listdir( options.data_dir ) if re.search( 'trj\d+\.npy', file )  ]
trjList.sort()
trjList = [ b for (a,b) in trjList ]

# The trj list is sorted now.
print "Loading Trajectories ..."
for trj in trjList:
#    print trj
    QList.extend( load( os.path.join( options.data_dir, trj ) ).astype(uint8)[::options.stride] )

QList = array( QList )
#QList = array([ bitarray(thing) for thing in QList ] )
print "K-Centers Clustering with %d conformations..." % QList.shape[0] 
a=time()
KCentersInd = Clustering.KCenters.ClusterQNorm( QList, options.numGens, QCutoff = options.cutoff, Seed= random.randint( len(QList ) ) )
b=time()
print b-a
save( 'Gens.npy', QList[ KCentersInd ] )
