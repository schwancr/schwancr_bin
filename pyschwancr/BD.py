#!/usr/bin/env python

import numpy as np
import random

class BrownianDynamicsIntegrator:
	""" This is a brownian dynamics integrator """
	def __init__(self,x0,ForceFunc,gammas=None,N=1,dt=1,noiseMean=0, noiseStdev=1, T=100, lower_bounds = None, upper_bounds = None):
		# Positions is initialized with the 
		self.positions = np.array([ x0 ]).astype(float)
		if gammas == None:
			self.gammas = np.ones( N )
		else:
			self.gammas = np.array([ gammas ])

		self.dt = dt

		self.dtOverGammas = float(dt) / self.gammas
		
		self.ForceFunc = ForceFunc

		self.particles = N

		self.noiseMean = noiseMean
		self.noiseStdev = noiseStdev

		self.temp = T
		self.randConstant = np.sqrt( 2. * 1.38E-23 * self.temp * self.dtOverGammas )
		if upper_bounds == None:
			upper_bounds = [ 1. ] * len(x0)
		if lower_bounds == None:
			lower_bounds = [ 0. ] * len(x0)

		self.upper_bounds = np.array( upper_bounds )
		self.lower_bounds = np.array( lower_bounds )
	def next(self):
		# Constants in the eqn:
		randG_x = random.normalvariate( self.noiseMean, self.noiseStdev )	
		randG_y = random.normalvariate( self.noiseMean, self.noiseStdev )
		randG = np.array( [ randG_x, randG_y ] )
		newShape = np.array( self.positions.shape )
		newShape[0] += 1 # Add one to the first index. So axis=0 corresponds to time.
		# But there is no restriction on how many other axes there are! In fact the entire structure 
		# is a black box since only F will need to understand it (which is user-defined)
		lastPositions = self.positions[-1].copy() # I make a copy because I need to modify this last entry
#		print lastPositions,
		self.positions.resize( newShape, refcheck = False )
		# Reshaped positions in order to add in my next iteration
		forces = self.ForceFunc( lastPositions )

		self.positions[-1] = lastPositions + self.dtOverGammas * forces + self.randConstant * randG

		#print '\t', self.dtOverGammas * forces + self.randConstant * randG, '\t', self.positions[-1]

		# Correct for going out of bounds:
		lowCheckArray = self.positions[-1] - self.lower_bounds
		upCheckArray = self.upper_bounds - self.positions[-1]
		#print self.positions[-1],
		self.positions[-1][ np.where( lowCheckArray < 0 ) ] = self.lower_bounds[ np.where( lowCheckArray < 0 ) ]
		self.positions[-1][ np.where( upCheckArray < 0 ) ] = self.upper_bounds[ np.where( upCheckArray < 0 ) ]
		#print '\t', self.positions[-1]

		#print self.dtOverGammas * forces, '\t', self.randConstant * randG
		#print '\t', self.dtOverGammas * forces + self.randConstant * randG
