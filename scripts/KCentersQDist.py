#!/home/schwancr/epd-7.1-1-rh5-x86_64/bin/python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='data_dir',help='Directory to look for numpy arrays of data (Q over time)')
parser.add_option('-u',type=int, default=1,dest='stride', help='Stride to sample by')
parser.add_option('-r',dest='cutoff',type=float, default=-1,help='Cutoff (% Similar NCs)')


options, args = parser.parse_args()

from numpy import *
from msmbuilder import Clustering
import os
import re

QList = []

trjList = [ ( int( file[3:-4] ), file ) for file in os.listdir( options.data_dir ) if re.search( 'trj\d+\.npy', file )  ]
trjList.sort()
trjList = [ b for (a,b) in trjList ]

# The trj list is sorted now.
print "Loading Trajectories ..."
for trj in trjList:
#    print trj
    QList.extend( load( os.path.join( options.data_dir, trj ) ).astype(int)[::options.stride] )

QList = array( QList )
print "K-Centers Clustering with %d conformations..." % QList.shape[0] 
KCentersInd = Clustering.KCenters.ClusterQDist( QList, 1000000, QCutoff = options.cutoff, Seed= random.randint( len(QList ) ) )

save( 'Gens.npy', QList[ KCentersInd ] )
