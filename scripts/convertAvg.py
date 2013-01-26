#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-f', dest='input_FN', default='./stateAvg_RMSD.Fixed.dat', help='Input stateAvg_RMSD.dat' )
options, args = parser.parse_args()
 
from msmbuilder import Serializer
from pyschwancr import dataIO
import os, sys, re
 
msmAvg = dataIO.readData( options.input_FN )

avgs = msmAvg[:,1]
vars = msmAvg[:,2] ** 2

s = Serializer.Serializer( { 'state_mean_rmsd': avgs, 'state_var_rmsd': vars } )

s.SaveToHDF( 'ClusterStats.hdf' )
