#!/usr/bin/env python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='data_FN',help='Data to fit (should be two columns, x and y)')
parser.add_option('-o',dest='out_FN',help='Output will write the average and standard deviations of the two gaussians')
parser.add_option('-b',dest='bins',type=int,default=100,help="Number of bins to use in the histogram")
parser.add_option('-w',dest='weights',help='Weights to use in creating the histogram')
options,args = parser.parse_args()

from numpy import *
from scipy import *
from scipy import optimize
from pyschwancr import dataIO
import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *


def fitFnc( p, x ):
	return p[0] * exp( - ( ( x - p[1] ) / p[2] )**2 ) + p[3] * exp( - ( ( x - p[4] ) / p[5] )**2 )

def errFnc( p, x, y):
   return (fitFnc( p, x ) - y) ** 2

def main():
	rawData = dataIO.readData(options.data_FN)
	if options.weights:
		weights = dataIO.readData( options.weights )
		sol = hist( rawData, bins = options.bins, weights = options.weights, label='Data' )
	else:
		sol = hist( rawData, bins = options.bins, label='Data' )

	xDat = array([ (sol[1][i]+sol[1][i+1])/2 for i in range( options.bins ) ])
	yDat = sol[0]
	data = array( zip( xDat, yDat ) )

	maxVal = data[:,1].max()
	minX = data[:,0].min()
	maxX = data[:,0].max()
	temp = (maxX - minX) * 0.1
	p0 = [ maxVal, temp, temp, maxVal, maxX - temp, temp ]
	
	pF, success = optimize.leastsq( errFnc, p0[:], args=( data[:,0], data[:,1] ), maxfev=100000, xtol=1E-10 )
#	print pF

	xFit = linspace(data[:,0].min(),data[:,0].max(),1000)
	yFit = fitFnc( abs(pF), xFit )

	#plot( data[:,0], data[:,1],'.',label='Data')
	plot( xFit, yFit, '-r',lw=2,label='Fit')
	legend()
	if options.out_FN:
		ttl = '.'.join( options.out_FN.split('.')[:-1] )
	else:
		ttl = '.'.join( options.data_FN.split('/')[-1].split('.')[:-1]) + 'Fit to Gaussians'

	xlabel('.'.join( options.data_FN.split('/')[-1].split('.')[:-1] ) )
	ylabel('Frequency')
	title(ttl)

	if options.out_FN:
		savefig(options.out_FN[:-4]+'.pdf')
	else:
		savefig('.'.join( options.data_FN.split('/')[-1].split('.')[:-1] ) + 'fit2gauss.pdf')

	print "Wrote parameters to %s" % dataIO.writeData( [ options.out_FN ], pF )

	return

if __name__ == '__main__':
	main()
