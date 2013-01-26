#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-V',dest='pot_FN',help="Energy surface to use for calculating gradients! FORMAT should be a 2D array with values for Log10(Counts). This will assume that the name looks like X_vs_Y.dat. Use plot2Dhist.py to generate this file")
parser.add_option('--xmin',dest='x_min',default=0.,type=float,help="Minumum in x direction")
parser.add_option('--ymin',dest='y_min',default=1.,type=float,help="Minumum in y direction")
parser.add_option('--xmax',dest='x_max',default=0.,type=float,help="Maximum in x direction")
parser.add_option('--ymax',dest='y_max',default=1.,type=float,help="Maximum in y direction")
parser.add_option('-g',dest='gamma',default=1.,type=float, help="Friction coefficient")
parser.add_option('-t',dest='temp',default=100.,type=float,help="Temperature")
parser.add_option('--dt',dest='dt',default=0.0005,type=float, help="Timestep")
parser.add_option('-n',dest='num_iter',default=100, type=int, help="Number of iterations to perform")
parser.add_option('-o',dest='out_FN',help="Output filename")
options, args = parser.parse_args()

from pyschwancr import BD, dataIO
from numpy import *
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
import random
import re

V = - 1.38E-23 * options.temp * np.log( 10**( dataIO.readData( options.pot_FN )  ) )
V = V.T # This is needed to swap the axes
Nx, Ny = V.shape
print V.shape
dx = 1. / (Nx - 1.)
dy = 1. / (Ny - 1.)

imshow(V.T,extent=[0,1,0,1],origin='bottom',cmap='jet_r') # V has x-axis = axis 0, but need to transpose to use imshow correctly

def Force( R ):
	wallForce = 100
	# Need to get the indices to look in V:
	I_x = int( R[0] * ( Nx - 1 ) )
	I_y = int( R[1] * ( Ny - 1 ) )

	if I_x >= Nx-1:
		I_x = Nx - 2
	if I_y >= Ny-1:
		I_y = Ny -2

	#if I_x < 0:
	#	print "BAD!", R, I_x, I_y
	#	I_x = 0
	#if I_y < 0:
	#	print "BAD!", R, I_x, I_y
	#	I_y = 0

	#if R[0] >= 1:
	#	F_x = - ( ( 1 - R[0] ) ** 2 + wallForce )
	#elif R[0] < 0:
	#	F_x = ( 0 - R[0] ) ** 2 + wallForce
	#else:
	F_x = - ( V[ I_x + 1 , I_y ] - V[ I_x, I_y ] ) / dx

	#if R[1] >= 1:
	#	F_y = - ( ( 1 - R[1] ) ** 2 + wallForce )
	#elif R[1] < 0:
	#	F_y = ( 0 - R[1] ) ** 2 + wallForce
	#else:
	F_y = - ( V[ I_x, I_y + 1 ] - V[ I_x, I_y ] ) / dy

	return array([ F_x, F_y ])
	

def main():
	

#	x0 = [ random.uniform(0,1), random.uniform(0,1) ]
	x0 = [ 0.1, 0.2 ]
	bd = BD.BrownianDynamicsIntegrator(x0, Force, dt = options.dt, T = options.temp, gammas=options.gamma, noiseStdev=1)

	for i in range(options.num_iter):
		bd.next()
		#print Force( bd.positions[-2] )
		#print bd.positions[-1]
	savetxt( 'test.dat', bd.positions )


	plot( bd.positions[:,0], bd.positions[:,1],color='white' )
	plot( x0[0], x0[1], '*',color='lightgreen',label="Start", markersize=10)
	plot( bd.positions[-1,0], bd.positions[-1,1], color='red',marker='o',ls='',label="End")
	colorbar().set_label( '-kT ln( Counts )' )
	#xlim([0,1])
	#ylim([0,1])
	legend( loc=2 )

	m = re.search( '(.*)_vs_(.*)\.dat', options.pot_FN )

	xlabel( m.group(1) )
	ylabel( m.group(2) )
	savefig(options.out_FN)



if __name__ == '__main__':
	main()


