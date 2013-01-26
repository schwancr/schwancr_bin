#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='dist_FN',default='./Assignments.h5.RMSD',help='Distance to generator from msmbuilder [ ./Assignments.h5.RMSD ]')
parser.add_option('-r',dest='cutoff',type=float, help='Cutoff used in clustering')
parser.add_option('--xl',dest='x_lbl',default='Distance to Generator',help='Label for the x-axis')
parser.add_option('-t',dest='title',default='Distance to Generator Distribution', help='Title for the plot')
parser.add_option('-o',dest='out_FN',default='Dist2Gen.pdf',help='Output filename to save the plot to')
options, args = parser.parse_args()
 
from numpy import *
from msmbuilder import Serializer
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from pyschwancr import dataIO
import os, sys, re

Dists = Serializer.LoadData( options.dist_FN )

Dists = Dists[ where( Dists >= 0 ) ]
np.savetxt( options.out_FN[:-4] + '.dat', Dists )
hist( Dists, bins=100 )

vlines( options.cutoff, 0, ylim()[1], color='red' )

xlabel( options.x_lbl )
ylabel( 'Frequency' )
title( options.title )
savefig( options.out_FN )
