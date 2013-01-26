#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()

parser.add_option('-U',dest="U_states",help="Unfolded states")
parser.add_option('-F',dest="F_states", help="Folded states")
parser.add_option('-t',dest="tProb",help="tProb.mtx from msmbuilder")
parser.add_option('-o',dest="out_FN",default="./FCs.dat",help="Output file to save forward committors")

options, args = parser.parse_args()

from scipy.io import mmread
from numpy import *
from msmbuilder import TransitionPathTheory as TPT

# Load data
Uind = loadtxt( options.U_states, int )
Find = loadtxt( options.F_states, int )
tProb = mmread( options.tProb )
print "Loaded data."
print "Calculating FCs"
Find = array( Find )
Uind = array( Uind )

if Find.shape == ():
	Find = array( [ Find ] )
if Uind.shape == ():
	Uind = array( [ Uind ] )
FC = TPT.GetFCommittors( Uind, Find, tProb, maxiter=10000 )
print "Done calculating, saving data to %s" % options.out_FN
savetxt( options.out_FN, FC )


