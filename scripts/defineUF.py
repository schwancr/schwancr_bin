#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-f',dest='data_FN',help='Average data for each state in your msm.')
parser.add_option('--low-is-folded',dest='low_is_folded',default=False,action='store_true',help='Pass this flag if low values in your order parameter indicate folded')
parser.add_option('--fc',dest='f_cut',type=float,help='Folded state cutoff')
parser.add_option('--uc',dest='u_cut',type=float,help='Unfolded state cutoff')
 
options, args = parser.parse_args()
 
from numpy import *
from msmbuilder import Serializer
from pyschwancr import dataIO
import os, sys, re
 
data = dataIO.readData( options.data_FN )
mObj = re.search('stateAvg_(.*).Fixed.dat', options.data_FN )
datName = mObj.group(1)

# First find the indices of the folded state:

avg_data = data[:,1]
std_data = data[:,2] 

if options.low_is_folded: # Then low means folded
	isFolded = lambda x : array([ i < options.f_cut for i in x ])
	isUnfolded = lambda x : array([ i > options.u_cut for i in x ])
else: # Then low means unfolded
	isFolded = lambda x : array([ i > options.f_cut for i in x ])
	isUnfolded = lambda x : array([ i < options.u_cut for i in x ])

avgPlusStd = avg_data + std_data
avgMinusStd = avg_data - std_data

print avgPlusStd.max(), avgPlusStd.min()
print avgMinusStd.max(), avgMinusStd.min()

FoldInd = where( isFolded( avgPlusStd ) * isFolded( avgMinusStd ) )[0]
UnfoldInd = where( isUnfolded( avgPlusStd ) * isUnfolded( avgMinusStd ) )[0]


print FoldInd, UnfoldInd

savetxt( 'F_states_%s.dat' % datName, FoldInd, '%d' )
savetxt( 'U_states_%s.dat' % datName, UnfoldInd, '%d' )
