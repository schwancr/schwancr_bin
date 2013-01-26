#!/home/schwancr/Installed/epd/bin/python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='data_FN',help='Data to fit (should be two columns, x and y)')
parser.add_option('-o',dest='out_FN',help='Output will write the average and standard deviations of the two gaussians')
parser.add_option('-w',dest='weight_FN',help='Weights if you want to add gaussians from stateAvg_<param>.dat')
options,args = parser.parse_args()

from numpy import *
from scipy import *
from scipy import optimize
import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *
from pyschwancr import dataIO


def fitFnc( p, x ):
   return p[0] * exp( - ( ( x - p[1] ) / p[2] )**2 ) + p[3] * exp( - ( ( x - p[4] ) / p[5] )**2 )

def errFnc( p, x, y):
   return (fitFnc( p, x ) - y) ** 2

def gauss( p, x ):
   return p[0] * exp( - ( ( x - p[1] ) / p[2] )**2 )

def prepData( data ):

   if len( data.shape ) == 1:
      return distribute( data )
   elif data.shape[1] == 1:
      return distribute( data[:,0] )
   elif data.shape[1] == 2:
      return data
   elif data.shape[1] == 3:
      return addGauss( data )
   else:
      print "Not sure how to handle this data. Shape should be one of: (N,); (N,1); (N,2); (N,3)"
      exit()

def addGauss( data ):
   if not options.weight_FN:
      print "Need weights in order to add gaussians together."
   
   means = data[:,1]
   stds = data[:,2]

   meanPlus2sigma = means + 2 * stds
   meanMinus2sigma = means - 2 * stds

   minX = min( [ 0, meanMinus2sigma.min() ] )
   maxX = meanPlus2sigma.max()
   print minX, maxX
   weights = dataIO.readData( options.weight_FN )

   weights *= sqrt( 2. * pi * stds )

   xi = linspace( minX, maxX, 100 )
   
   y = zeros( 100 )

   for i in range( data.shape[0] ):
      y += gauss( [ weights[i], means[i], stds[i] ], xi )
   
   return array( zip( xi, y ) )

def distribute( data ):
   y, x0, c = hist( data, bins = 100 )
   x = [ ( x0[i] + x0[i+1] ) / 2. for i in range( len( x0 ) - 1 ) ]

   return array( zip(x,y) )

def main():
	data = dataIO.readData(options.data_FN)
	data = prepData( data )
	maxVal = data[:,1].max()
	minX = data[:,0].min()
	maxX = data[:,0].max()
	temp = (maxX - minX) * 0.1
	p0 = [ maxVal, temp, temp, maxVal, maxX - temp, temp ]
	
	pF, success = optimize.leastsq( errFnc, p0[:], args=( data[:,0], data[:,1] ), maxfev=100000, xtol=1E-10 )
	print pF

	xFit = linspace(data[:,0].min(),data[:,0].max(),100)
	yFit = fitFnc( abs(pF), xFit )

	plot( data[:,0], data[:,1],'.',label='Data')
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
