#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='data_FN',help='Data to define F and U ensembles. Should be an Nx3 matrix corresponding to index, average, standard deviation. ')
parser.add_option('-n',dest='num',type=int,default=50,help='Number of states to sample from F and U')
parser.add_option('-c',dest='cut',type=float,help='Cutoff to use to define F and U ensembles based on the data.')
parser.add_option('-p',dest='pops_FN',default='./clusterPops.dat',help='Cluster populations (use getClusterPops.py) [ ./clusterPops.dat ]')
parser.add_option('-m',dest='minCount',type=int,default=10,help='Minimum counts in the state to use it [ 10 ]')
parser.add_option('-o',dest='out_FN',default='States.dat',help='Output filename [ States.dat ]')
parser.add_option('--nf',dest='num_out',default=1,type=int,help='Number of output files to write')
parser.add_option('--add-most-populated',dest='add_most_pop',default=False,action='store_true',help='Add the most populated state to the list, and remove one from that ensemble')
options, args = parser.parse_args()
 
import numpy as np
from pyschwancr import dataIO
import os, sys, re
 
print "Loading data..."

Data = dataIO.readData( options.data_FN )[:,1]
Pops = dataIO.readData( options.pops_FN )
outFN = options.out_FN
cut = options.cut
num = options.num
minCount = options.minCount

print "Working..."

MostPopInd = np.argmax( Pops )
MostPopIsAbove = ( Data[ MostPopInd ] > cut )

EnoughInd = np.where( Pops >= minCount )[0]
EnoughData = Data[ EnoughInd ]

BelowStates = EnoughInd[ np.where( EnoughData < cut ) ]
AboveStates = EnoughInd[ np.where( EnoughData > cut ) ]

RandBelow = np.random.permutation( BelowStates )[:num]
RandAbove = np.random.permutation( AboveStates )[:num]

if options.add_most_pop: # If we want to add the most populated state, then
    if MostPopIsAbove: # If the most populated state is above the cut, then
        if not MostPopInd in RandAbove: # If the most pop state is not already in the state list
            RandAbove[-1] = MostPopInd # replace the last state with the most pop state
    else: # If the most populated state is below the cut
        if not MostPopInd in RandBelow: # If the most pop state is not already in the state list
            RandBelow[-1] = MostPopInd # replace the last state with the most pop state


States = np.concatenate( ( RandBelow, RandAbove ) )
States.sort()

statesPer =  int( States.shape[0] / float( options.num_out ) + 1 )
outFN = options.out_FN

for i in range( options.num_out ):
   if outFN[-4:] in [ '.dat','.txt','.npy']:
      outFN = outFN[:-4]
   
   np.savetxt( outFN + str(i) + '.dat', States[ i * statesPer : ( i + 1 ) * statesPer ], "%d")

