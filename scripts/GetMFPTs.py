#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-P',dest='pops_FN',default='./Populations.dat',help='Populations from msmbuilder [ ./Populations.dat ]')
parser.add_option('-T',dest='T_FN',default='./tProb.mtx',help='Transition probability matrix from msmbuilder [ ./tProb.mtx ]')
parser.add_option('-l',dest='lag',type=float,help='Lag time used to build T')
parser.add_option('-o',dest='out_FN',default='MFPTs.npy',help='MFPTs for states specified')
parser.add_option('-S',dest='state_FN',help='File with a list of states to calculate the MFPT TO. If no file is specified, then ALL states will be used')
parser.add_option('-p',dest='procs',default=1,type=int,help='Number of processes to run [ 1 ]')
options, args = parser.parse_args()
 
import numpy as np
from scipy.io import mmread
from msmbuilder import TransitionPathTheory as TPT
from pyschwancr import dataIO
import os, sys, re
 
def AnalyzeIndex( state ):
    print "Working on state %d" % state
    return TPT.GetMFPTsolve( [ state ], T, LagTime=Lag )
    #return TPT.GetMFPTFundMat( state, T, Pops, LagTime=Lag )

print "Loading data..."
T = mmread( options.T_FN )
Pops = dataIO.readData( options.pops_FN )
Lag = options.lag
outFN = options.out_FN
if options.state_FN:
    F = dataIO.readData( options.state_FN ) # If there is a list of states use it
else: # Otherwise check the MFPTs for ALL states in T
    F = np.arange( T.shape[0] )

print "Calculating MFPTs"

sol = []
for state in F:
   sol.append( AnalyzeIndex( state ) )

ResAry = np.array( sol ).T
# This is the results array. The data will be stored as columns, so the i,j th entry will be the MFPT from state i to state j
print "Saving data to %s. Remember that this is a matrix who's columns are the MFPTs to the corresponding state in the input state set." % outFN

np.save( outFN, ResAry )
