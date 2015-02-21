import numpy as np

class LangevinDynamicsIntegrator():

    def __init__( self, V, dV, dt = 0.001, gamma = 1, kT = 0.25, xRange=None ):
    	self.V = V
    	self.dV = dV
    	self.gamma = gamma
    	self.kT = kT
    	self.positions = np.array([])
    	self.randCnst = np.sqrt( 2 * kT * gamma )
    	self.dt = dt 
    	self.xRange = np.array( xRange )
    	self.stdev = np.sqrt( self.dt )
    	return
    
    def start( self, x0 ):
    	self.positions = np.array([ x0 ])	
    	return

    def next( self ):
    	randG = np.random.normal(0,self.stdev, size=2) # Default is mean = 0 stdev = 1

    	xi = self.positions[-1]

        #print "xi:", xi
    	xf = xi + ( - self.dV( xi ) * self.dt + self.randCnst * randG ) / self.gamma
    	#print "xf:", xf

    	too_small = np.where( xf < self.xRange[0] )
    	xf[ too_small ] = 2 * self.xRange[0][ too_small ] - xf[ too_small ]
        #print "small:", xf

    	too_big = np.where( xf > self.xRange[1] )
    	xf[ too_big ] = 2 * self.xRange[1][ too_big ] - xf[ too_big ]
        #print "big:", xf
    	#if xf < self.xRange[0]:
    #		xf = self.xRange[0] + ( self.xRange[0] - xf )
    #	elif xf > self.xRange[1]:
    #		xf = self.xRange[1] + ( self.xRange[1] - xf )

    	self.positions = np.concatenate( (self.positions, [ xf ] ) )

    	return
    
