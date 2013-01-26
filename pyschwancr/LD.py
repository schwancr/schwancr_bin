import numpy as np

class LangevinDynamicsIntegrator():

	def __init__( self, V, dV, dt = 0.001, gamma = 1, kT = 0.25, xRange = np.array([0,1]) ):
		self.V = V
		self.dV = dV
		self.gamma = gamma
		self.kT = kT
		self.positions = np.array([])
		self.randCnst = np.sqrt( 2 * kT * gamma )
		self.dt = dt 
		self.xRange = xRange
		self.stdev = np.sqrt( self.dt )
		return
	
	def start( self, x0 ):
		self.positions = np.array([ x0 ])	
		return

	def next( self ):
		randG = np.random.normal(0,self.stdev) # Default is mean = 0 stdev = 1

		xi = self.positions[-1]
		xf = xi + ( - self.dV( xi ) * self.dt + self.randCnst * randG ) / self.gamma
		
		if xf < self.xRange[0]:
			xf = self.xRange[0] + ( self.xRange[0] - xf )
		elif xf > self.xRange[1]:
			xf = self.xRange[1] + ( self.xRange[1] - xf )

		self.positions = np.concatenate( (self.positions, [ xf ] ) )

		return
	
