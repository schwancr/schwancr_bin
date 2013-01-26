"""
This module will do montecarlo moves for a 2 dimensional lattice.
"""
import numpy as np


def checkDiag( mover, rTest ):

		if len( mover.positions.shape ) < 2:
			start = mover.positions
		else:
			start = mover.positions[-1]
		startCount = mover.Pot[ start[0], start[1] ]
		testCount = mover.Pot[ tuple( rTest ) ]

		# Need to test if the counts are zero.
		if startCount == testCount:
			return True
		elif startCount == 0:
			return True
		

		prob = min([ 1, testCount / float( startCount ) ] )

		if np.random.random() < prob:
			return True
		else:
			return False

class MC_Mover:
	def __init__(self, xData, yData, xRange = [0,1], yRange = [0,1], takeMoveFunc = checkDiag ):
		"""
		This function initializes the montecarlo mover
		
		"""
		# First get the number of levels on each axis.
		self.Nx = 1. / ( np.unique(xData)[1:] - np.unique(xData)[:-1] ).min() + 1
		self.Ny = 1. / ( np.unique(yData)[1:] - np.unique(yData)[:-1] ).min() + 1
		# Adding one is necessary to include both boundaries

		self.Pot = np.zeros( ( self.Nx, self.Ny ) ) 
		# This means the ROWS are the x-axis and the COLUMNS are the y-axis
		# So for imshow, it will need to be transposed to look right.

		xPts = xData * (self.Nx - 1) / ( xRange[1] - xRange[0] )
		yPts = yData * (self.Ny - 1) / ( yRange[1] - yRange[0] )
		self.Nx = int( self.Nx )
		self.Ny = int( self.Ny )
		for i in xrange( len( xData ) ):
			self.Pot[ xPts[i], yPts[i] ] += 1

		# The above defines the counts potential
		self.xTicks = np.linspace( xRange[0], xRange[1], self.Nx )
		self.yTicks = np.linspace( yRange[0], yRange[1], self.Ny )
		
		self.positions = np.array([])
		
		self.takeMove = takeMoveFunc

		return

	def start( self, r0 ):
		# This script starts the iteration at some beginning value. (It can also restart the simulation by deleting the old positions.)
		self.positions = np.array( r0 )
		self.positions.resize( ( 1, len( r0 ) ) )
		return

	def next( self ):
		
		if self.positions == np.array([]):
			print "ERROR: Need to begin the simulation with self.start( r0 )."
			exit()
		
		moveAry = np.array([ -1 if i==0 else 1 for i in np.random.randint(0,2,2)] )
		# The move array will be a 0 for decrease and 1 for increase
		rTest = self.positions[-1] + moveAry

		rTest = np.array([ 0 if coord < 0 else coord for coord in rTest ])
		# Above line makes sure neither coordinate is below 0.
		rTest = np.array([ (self.Nx-1, self.Ny-1)[i] if rTest[i] > (self.Nx-1, self.Ny-1)[i] else rTest[i] for i in (0,1) ] )
		# Above line makes sure neither coordinate is above the number of entries on that axis.

		if self.takeMove( self, rTest ):
			self.positions = np.vstack( ( self.positions, rTest ) )
		else:
			self.positions = np.vstack( ( self.positions, self.positions[-1] ) )

		return

	def getXYindices( self, r):
		
		# Find the closest thing in self.xTicks
		
		minXary = abs( self.xTicks - r[0] )
		minYary = abs( self.yTicks - r[1] )
		
		xInd = np.where( minXary == minXary.min() )[0][0]
		yInd = np.where( minYary == minYary.min() )[0][0]

		return np.array([ xInd, yInd ])

	def getXYvalues( self ):
		xDat = self.xTicks[ self.positions[:,0] ]
		yDat = self.yTicks[ self.positions[:,1] ]

		return np.array( zip( xDat, yDat ) )
		
def RandomWalk( MC_mover, r0, N=100, untilFcn=None):
	"""
	This function will do a random walk on a montecarlo mover object given an initial state and number of iterations
	"""

	if untilFcn == None:
		untilFcn = lambda r : False
		# i.e. Never stop

	MC_mover.start( r0 )
	
	for i in xrange( N ):
		MC_mover.next()
		if untilFcn( MC_mover.positions[-1] ):
			break

	return untilFcn( MC_mover.positions[-1] )	
	


