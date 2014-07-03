#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-a','--assignment-file',dest='assFN',default='./Assignments.h5',help='Assignments.h5 file generated by msmbuilder')
parser.add_option('-o','--output-file',dest='outFN',default='./clusterPops.dat',help='Flat text file with the population of each cluster')

options, args = parser.parse_args()

import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *
from msmbuilder import io
from numpy import *

# First read in the data:

As = io.loadh(options.assFN)
try:
    As = As['Data']
except:
    As = As['arr_0']

As1D = As[ where( As >= 0 ) ].flatten()

stuff = bincount( As1D )

savetxt(options.outFN,stuff, "%d")

b = bincount(stuff)

plot(b)
xscale('log')
ylabel('Frequency')
xlabel('Number of Assignments per State')

savefig( options.outFN[:-4]+'.pdf')
