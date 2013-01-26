#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option('-r',dest='raw_FN',help='Raw data for all the trajectories to calculate the autocorrelations')
parser.add_option('-o',dest='out_FN',help='Output filename to write to. Name used for the data and the plots')
parser.add_option('-P',dest='procs',default=1,type=int,help='Number of processes to use [ 1 ]' )
parser.add_option('--n-rates',dest='nRates',default=100,type=int,help='Number of rates to use in calculating the rate spectra [ 100 ]')
parser.add_option('--plot',dest='numPlot',type=int,default=3,help='Number of autocorrelations and their fits to plot')
parser.add_option('--times',dest='times',default='10,1e5',help='Timescales to look at. Format: --times a,b [ 10,1e5 ]')

options, args = parser.parse_args()
 
from msmbuilder import Serializer, autocorrelate
from pyschwancr import dataIO, msmTools
import os, sys, re
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import warnings 
from ratespec import RateSpecClass as RS
from ratespec import RateSpecTools as RStools
import multiprocessing as mp

#warnings.filterwarnings('ignore','Warning: overflow encountered in exp')
#warnings.filterwarnings('ignore','divide\sby\szero')

np.seterr( divide='ignore' )

def f( x , a, b ):
	return np.exp( - b * x )

def AnalyzeTraj( auto ):

	data = np.array( zip( np.arange( len(auto) ) , auto ) )
	rs = RS.RateSpecClass( data = data, timeUnit = 1, nRates = numRates, RateRange=RateRange, standardizeData=False, scaleData=False, Lnorm='lasso' )
	A,rss,ndof = RStools.fitRateSpectrum( rs.Times, rs.Data, rs.Rates, 1.0E-3, Lnorm = rs.Lnorm, standardizeData=False )

	return A

numRates = options.nRates

print "Loading data"
proj = Serializer.Serializer.LoadFromHDF( options.proj_FN )

data = dataIO.readData( options.raw_FN )
data2d = msmTools.reshapeRawData( data, proj )

print "Calculating autocorrelations"
Autos = [ autocorrelate.fft_autocorrelate( trj[ np.where( trj != -1 ) ] ) for trj in data2d ]

print "Fitting the data to single exponentials"


TimeScales = []

RateRange = [ 1 / float( time ) for time in options.times.split(',') ]
RateRange.reverse()

pool = mp.Pool( options.procs )
result = pool.map_async( AnalyzeTraj, Autos )
result.wait()
sol = result.get()

# Since all of the autocorrelations are done in the same way rs.Timescales are the same for all of them. so just make a new one and get the timescales.
rs = RS.RateSpecClass( data = np.array( zip( np.arange( len( Autos[0] ) ), Autos[0] ) ), timeUnit = 1, nRates = numRates, RateRange=RateRange, standardizeData=False, scaleData=False, Lnorm='lasso' )

Spectra = np.array( sol )
TimeScales = rs.Timescales[ np.argmax( Spectra, axis = 1 ) ]

TimeScaleRange = rs.Timescales

#for i, auto in enumerate( Autos ):
#	print "\rWorking on trajectory %d ..." % i,
#	data = np.array( zip( np.arange( len(auto) ) , auto ) )
#	rs = RS.RateSpecClass( data = data, timeUnit = 1, nRates = 2, RateRange=RateRange, standardizeData=False, scaleData=False, Lnorm='lasso' )
#	A,rss,ndof = RStools.fitRateSpectrum( rs.Times, rs.Data, rs.Rates, 1.0E-3, Lnorm = rs.Lnorm, standardizeData=False )

#	timeScale = np.argmax( A )
#	TimeScales.append( timeScale )

#	Spectra[i] = A

outName = '.'.join( options.out_FN.split('.')[:-1] )
TimeScales = np.array( TimeScales )
np.savetxt( outName + 'Timescale.dat', TimeScales )

Autos1D = np.concatenate( Autos )
np.save( outName + '.Autos.npy', Autos1D )

np.save( outName + '.RateSpec.npy', Spectra )
np.savetxt( outName + '.RateSpec_x.npy', TimeScaleRange )

print "Plotting some fits"

indList = np.random.permutation( np.arange( len( Autos ) ) )

minLength = min( [ len( trj ) for trj in Autos ] )

xi = np.linspace(0,minLength,10000)

Nplots = min( [ options.numPlot, 10 ] )

Colors = [ 'red', 'blue','purple','yellow','black','gray','cyan','navy','pink','green' ]
Markers = [ '.','o','s','^','*','+','x','D','<','>' ]
for i in range( Nplots ):
	ind = indList[i]

	name = 'Traj %d' % ind
	yi = f( xi, 1, 1 / TimeScales[i] )
	corr = Autos[ind]
	corr = corr[ np.where( corr != -1 ) ]

	plt.plot( range( len( corr ) ), corr, color = Colors[i],linestyle='-', label=name )
	plt.plot( xi, yi, label= name + ' Fit', color = Colors[i], linestyle='--' )

plt.xscale('log')
plt.ylim([0,1])
plt.legend()
plt.savefig( outName + '.png' )
