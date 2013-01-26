import numpy as np
from scipy import optimize

np.seterr( over='ignore' )
# Functions to use for the fitting process
def gauss( p, x):
	sol = np.sqrt( p[0] ** 2 ) * np.exp( - 0.5 * ( (x-p[1])/p[2] )**2 ) 
	return sol

def multiGauss( p, x ):
	pNew = p.reshape( ( p.shape[0]/3, 3 ) )
	sol = 0.
	for params in pNew:
		sol += gauss( params, x )
	return sol
	
def errFuncMultiGauss( p, x, y ):
	sol = multiGauss(	p, x ) - y 
	return sol

def GaussFit( xDat, yDat, N=1, PadWithZeros = True ):
	""" 
This function fits data to a sum of gaussians. 

Input:
	1) xDat - np.array with x values
	2) yDat - np.array with y values to fit to
	3) N [ 1 ] - Number of gaussians to fit (as a sum)

Output:
	1) p - np.array with the parameters. Each triplet corresponds to that particular gaussians height, mean, standard deviation
	"""


	#p0 = np.ones( N * 3 )
	p0 = np.array([ yDat.max(), xDat.mean(), xDat.std() ] )
	# This will determine how many parameters there are.
	xRange = xDat.max() - xDat.min()
	Npts = len( xDat )
	if PadWithZeros:
		x0 = xDat.min() - 10 * xRange
		x1 = xDat.min() - xRange / Npts
		x2 = xDat.max() + xRange / Npts
		x3 = xDat.max() + 10 * xRange

		X = np.concatenate( ( np.linspace( x0, x1, Npts-1 ), xDat, np.linspace( x2, x3, Npts-1 ) ) )
		Y = np.concatenate( ( np.zeros( Npts-1 ), yDat, np.zeros( Npts-1 ) ) )
		# Pad the data with zeros outside the range of the data
	else:
		X = xDat
		Y = yDat

	sol = optimize.leastsq( errFuncMultiGauss, x0 = p0, args = ( X, Y ), maxfev = 10000000 )
	if sol[1] in [1,2,3,4]:
		return sol[0]
	else:
		print "Error in least squares fit"
		return -1

# THE FOLLOWING ARE GAUSSIANS PLUS A CONSTANT
# Functions to use for the fitting process 

def gaussPlusC( p, x):
   sol = np.sqrt( p[0] ** 2 ) * np.exp( - ( (x-p[1])/p[2] )**2 ) + p[3] 
   return sol

def multiGaussPlusC( p, x ):
   pNew = p.reshape( ( p.shape[0]/4, 4 ) )
   sol = 0.
   for params in pNew:
      sol += gaussPlusC( params, x )
   return sol

def errFuncMultiGaussPlusC( p, x, y ):
   sol = multiGaussPlusC( p, x ) - y
   return sol

def GaussFitPlusC( xDat, yDat, N=1, PadWithZeros = True ):
   """ 
This function fits data to a sum of gaussians. 

Input:
   1) xDat - np.array with x values
   2) yDat - np.array with y values to fit to
   3) N [ 1 ] - Number of gaussians to fit (as a sum)

Output:
   1) p - np.array with the parameters. Each triplet corresponds to that particular gaussians height, mean, standard deviation
   """


   #p0 = np.ones( N * 3 )
   p0 = np.array([ yDat.max(), xDat.mean(), xDat.std(), yDat.min() ] )
   # This will determine how many parameters there are.
   xRange = xDat.max() - xDat.min()
   Npts = len( xDat )
   if PadWithZeros:
      x0 = xDat.min() - 10 * xRange
      x1 = xDat.min() - xRange / Npts
      x2 = xDat.max() + xRange / Npts
      x3 = xDat.max() + 10 * xRange

      X = np.concatenate( ( np.linspace( x0, x1, Npts-1 ), xDat, np.linspace( x2, x3, Npts-1 ) ) )
      Y = np.concatenate( ( np.zeros( Npts-1 ), yDat, np.zeros( Npts-1 ) ) )
      # Pad the data with zeros outside the range of the data
   else:
      X = xDat
      Y = yDat

   sol = optimize.leastsq( errFuncMultiGaussPlusC, x0 = p0, args = ( X, Y ), maxfev = 10000000 )
   if sol[1] in [1,2,3,4]:
      return sol[0]
   else:
      print "Error in least squares fit"
      return -1
# Exponential functions

def expon( p, x):
	sol = p[0] * np.exp( - x / p[1] )
	return sol

def multi_expon( N, p, x ):
    const = p[-1]
    p = p[:-1].reshape( (N,-1) )
    
    sum = const
    for i in range( N ):
        sum += expon( p[i], x )

    return sum
    
def exponHt1( p, x ):
	sol = np.exp( - p[0] * x )
	return sol

def errFuncMultiExpon( p, x, y, N ):
	sol = multi_expon( N, p, x) - y
	return sol

def errFuncExpon( p, x, y ):
	sol = expon( p, x) - y
	return sol

def errFuncExponHt1( p, x, y ):
	sol = exponHt1( p, x) - y
	return sol

def ExponFit( xDat, yDat, LogSample = False, nPts = 1000 ):
	"""
This function fits data to an exponential. It does this with a least squares regression,
	however you can first subsample the data with a log scale by passing LogSample = True.

Inputs:
	1) xDat - np.array with x values
	2) yDat - np.array with y values
	3) LogSample [ False ] - Pass this if you want to subsample the data to remove some points at large x.

Outputs:
	1) p - np.array with the parameters to this fit: p[0] * np.exp( - p[1] * x )
	"""

	if LogSample:
		ind = np.logspace( 0, np.log10( xDat.shape[0] - 1 ), nPts )
		ind = np.unique( ind ).astype(int)
		xDat = xDat[ ind ]
		yDat = yDat[ ind ]

	sol = optimize.leastsq( errFuncExponHt1, x0 = [0.1], args = (xDat, yDat ),maxfev = 1000000 )
	if sol[1] in [1,2,3,4]:
		return sol[0]
	else:
		print "Error in least squares fit"
		return -1

def MultiExponFit( xDat, yDat, num_expons, LogSample = False, nPts = 1000, start_x0 = None ):
    """
This function fits data to an exponential. It does this with a least squares regression,
	however you can first subsample the data with a log scale by passing LogSample = True.

Inputs:
	1) xDat - np.array with x values
	2) yDat - np.array with y values
	3) LogSample [ False ] - Pass this if you want to subsample the data to remove some points at large x.

Outputs:
	1) p - np.array with the parameters to this fit: p[0] * np.exp( - p[1] * x )
    """

    if LogSample:
        ind = np.logspace( 0, np.log10( xDat.shape[0] - 1 ), nPts )
        ind = np.unique( ind ).astype(int)
        xDat = xDat[ ind ]
        yDat = yDat[ ind ]

    if start_x0 == None:
        start_x0 = np.ones( 2 * num_expons + 1 )
    
    sol = optimize.leastsq( errFuncMultiExpon, start_x0, args = ( xDat, yDat, num_expons ),maxfev = 1000000, xtol=1E-13, ftol=1E-13, full_output=True)
    if sol[-1] in [1,2,3,4]:
        print sol[-2]
        return sol[0]
    else:
        print "Error in least squares fit"
        return -1

