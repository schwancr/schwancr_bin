#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-c',dest='clust_FN', default='./clusterPops.dat',help='Cluster populations file from getClusterPops.py [ ./clusterPops.dat ]' )
parser.add_option('-o',dest='out_FN',default='./clusterPops.pdf',help='Output filename [ ./clusterPops.pdf ]')

options,args = parser.parse_args()

from numpy import *
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
import os
import re

# Load the data.

pops = loadtxt( options.clust_FN, int )

tots = bincount( pops )

plot( tots, '.' )
xscale( 'log' )
yscale( 'log' )
xlabel( 'Cluster population' )
ylabel( 'Frequency (# of clusters)' )

vlines( pops.mean(), ylim()[0], ylim()[1], color='red', label='Average Population' )


### Get plot title from the path name
pathName = os.path.abspath('.')
m = re.search( '/(Data.*)', pathName)
if m:
	title( "Cluster Population Distribution ( %s )" % m.group(1) )
else:
	title( "Cluster Population Distribution" )

text( 10 ** ( 0.7 * log10( xlim()[1] ) ), 10 ** ( 0.9 * log10( ylim()[1] ) ), "%s = %.1f\n%s = %d" % ( "Mean" , pops.mean(), "N", pops.shape[0] ) )
savefig( options.out_FN )

