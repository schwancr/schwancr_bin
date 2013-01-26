#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-x',dest='x_data',help='x-axis data')
parser.add_option('-y',dest='y_data',help='y-axis data')

options,args = parser.parse_args()
import numpy as np
from pyschwancr import dataIO
from scipy import optimize

def f( p, x ):
	return p[0] * x + p[1]

def err_f( p, x, y ):
	return f( p, x ) - y 

xDat = dataIO.readData( options.x_data )
yDat = dataIO.readData( options.y_data )

if len( xDat.shape ) > 1 :
	xDat = xDat[:,0]
	print "Using first column of x-axis data"

if len( yDat.shape ) > 1 :
	yDat = yDat[:,0]
	print "Using first column of x-axis data"

xMax = xDat.max()
xMin = xDat.min()
#xDat = ( xDat - xMin ) / xMax

yMax = yDat.max()
yMin = yDat.min()
#yDat = ( yDat - yMin ) / yMax

p0 = [1,1]
print "X-shape: %s, Y-shape: %s" % (str( xDat.shape ), str( yDat.shape ) )
sol = optimize.leastsq( err_f, x0 = p0, args = ( xDat, yDat ), ftol=1000)

m = sol[0][0]
b = sol[0][1]

print "m=%.4e" % m
print "b=%.4e" % b

np.savetxt('FitData_Y-%s_X-%s.dat' % (options.y_data[:-4].split('/')[-1] , options.x_data[:-4].split('/')[-1] ) , sol[0] )
