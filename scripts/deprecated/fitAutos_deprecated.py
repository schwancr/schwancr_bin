#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option('-r',dest='raw_FN',help='Raw data for all the trajectories to calculate the autocorrelations')
parser.add_option('-o',dest='out_FN',help='Output filename to write to. Name used for the data and the plots')
parser.add_option('--plot',dest='numPlot',type=int,default=3,help='Number of autocorrelations and their fits to plot')

options, args = parser.parse_args()
 
from msmbuilder import Serializer, autocorrelate
from pyschwancr import dataIO, msmTools
import os, sys, re
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import warnings 


warnings.filterwarnings('ignore','Warning: overflow encountered in exp')

def f( x , a, b ):
	return np.exp( - b * x )
print "Loading data"
proj = Serializer.Serializer.LoadFromHDF( options.proj_FN )

data = dataIO.readData( options.raw_FN )
data2d = msmTools.reshapeRawData( data, proj )

print "Calculating autocorrelations"
Autos = [ autocorrelate.fft_autocorrelate( trj[ np.where( trj != -1 ) ] ) for trj in data2d ]

print "Fitting the data to single exponentials"
Fits = [ curve_fit( f, np.arange( len( corr ) ), corr )[0] for corr in Autos ]

outName = '.'.join( options.out_FN.split('.')[:-1] )
Fits = np.array( Fits )
np.savetxt( outName + '.dat', Fits )


print "Plotting some fits"

indList = np.random.permutation( np.arange( len( Autos ) ) )

minLength = min( [ len( trj[ np.where( trj != -1 ) ] ) for trj in Autos ] )
xi = np.linspace(0,minLength,10000)

Nplots = min( [ options.numPlot, 10 ] )

Colors = [ 'red', 'blue','purple','yellow','black','gray','cyan','navy','pink','green' ]
Markers = [ '.','o','s','^','*','+','x','D','<','>' ]
for i in range( Nplots ):
	ind = indList[i]

	name = 'Traj %d' % ind
	yi = f( xi, 1, Fits[ind,1] )
	corr = Autos[ind]
	corr = corr[ np.where( corr != -1 ) ]

	plt.plot( range( len( corr ) ), corr, color = Colors[i],linestyle='-', label=name )
	plt.plot( xi, yi, label= name + ' Fit', color = Colors[i], linestyle='--' )

plt.xscale('log')
plt.ylim([0,1])
plt.legend()
plt.savefig( outName + '.png' )
