from numpy import *
import numpy as np
import multiprocessing as mp
from time import time 
from msmbuilder import Trajectory 
from pyschwancr import FitData

def calcPhi( U, T, F, NCs, pops):
   """
This function calculates the phi value for a particular residue (given by NCs).

Input:
   1) U - np.array of indices for the unfolded state
   2) T - np.array of indices for the TSE
   3) F - np.array of indices for the folded state
   4) NCs - np.array of averages for each state of NCs of a single residue
   5) pops - np.array of the populations of ALL states
Output:
   1) phi - float containing the phi value
"""

   phi = dot( NCs[T], pops[T] ) - dot( NCs[U], pops[U] )
   phi = phi / ( dot( NCs[F], pops[F] ) - dot( NCs[U], pops[U] ) )

   return phi

def calcFracFold( F, tProb, X0, N=10000 ):
   """
This function will calculate the fraction folded of an MSM over some number of timesteps.

Input:
   1) F - np.array of indices for the folded state
   2) tProb - scipy.sparse matrix in CSR format
   3) X0 - Initial populations vector
   4) N - Number of itrations to do

Outpout:
   1) fracFold - np.array of The fraction folded over time
   """
   sumFold = 0.
   X = X0
   fracFold = []
   for i in xrange(N):
      sumFold += X[F].sum()
      fracFold.append( sumFold )
      X[F] = 0.
      X = X * tProb
   return array( fracFold )

def getDegrees( T ):
	"""
This function gets the degrees of all states in T.
	
Input:
	1) T - scipy.sparse matrix corresponding to the *symmetric transition probability matrix

Output:
	1) Deg - np.array of degrees for each ROW in T
	"""

	rows, cols = T.nonzero()

	rowDeg = bincount( rows )
	colDeg = bincount( cols )

	if (rowDeg != colDeg).all():
		print "Row degree not the same as the column degree... This means your matrix is NOT symmetric"
		print "Returning the row degree..."


	return rowDeg

def calcRawFoldTime( traj, Fcut, Ucut, low_is_folded = True ):
	"""
This function will calculate the folding and unfolding time for a raw trajectory

Input:
	1) traj - np.array of data values as a function of time
	2) Fcut - cutoff defining the folded state
	3) Ucut - cutoff defining the unfolded state
	4) low_is_folded [ True ] - boolean specifying whether low or high is folded. Default is True (i.e. RMSD)

Output:
	1) FoldTime - np.array of folding times for all transitions in traj
	2) UnfoldTime - np.array of unfolding times for all transitions in traj
	"""
	if len( traj.shape ) != 1:
		print "Only one trajectory at a time please..."
		exit()
	
	if low_is_folded:
		isFolded = lambda x : x <= Fcut
		isUnfolded = lambda x : x >= Ucut
	else:
		isFolded = lambda x : x >= Fcut
		isUnfolded = lambda x : x <= Ucut

	states = ones( len( traj ) ) * -1 
	states[ where( isFolded( traj ) ) ] = 1
	states[ where( isUnfolded( traj ) ) ] = 0

	lastAss = -1
	for i in range( len( traj ) ):
		if states[i] == -1:
			states[i] = lastAss
		else:
			lastAss = states[i]

	# Now count the time spent in F or U

	transitions = logical_xor( states[:-1], states[1:] ) # This will be True at a transition!

	whereTrans = concatenate( ( [0], where( transitions )[0] ) )
	
	diffs = array([ whereTrans[i+1] - whereTrans[i] for i in range( len( whereTrans ) - 1 ) ])

	firstAss = min( where( states != -1 )[0] )

	if isFolded( states[ firstAss ] ):
		Folds = diffs[ firstAss : : 2 ]
		Unfolds = diffs[ firstAss+1 : : 2 ]
	else:
		Folds = diffs[ firstAss+1 : : 2 ]
		Unfolds = diffs[ firstAss : : 2 ]

	return Folds, Unfolds

def calcRg( traj ):
	"""
This function calculates the radius of gyration for a trajectory (Trajectory.Trajectory object)

Input
	1) traj - msmbuilder.Trajectory.Trajectory object

Output
	1) Rgs - np.array of Rgs over time for the trajectory
	"""

	Rgs = []
	traj.RestrictAtomIndices( where( traj['AtomNames'] == 'CA' )[0] )
	for frame in traj['XYZList']:
		rMean = frame.mean( axis = 0 )

		diffs = frame - rMean
		Rgs.append( ( diffs * diffs ).sum() )

	Rgs = array( Rgs )
	Rgs = Rgs / traj['XYZList'].shape[1]
	Rgs = sqrt( Rgs )
		
	
	return Rgs

def reshapeRawData( data, proj ):
	"""
	This function reshapes a 1d array of data corresponding to the data of each conformation. The shape is the same as assignments.h5

	Input:
	1) data - np.array of data to reshape
	2) proj - Project.Project object or dictionary with the key 'TrajLengths'

	Output:
	1) data2d - np.array of the correct shape.

	"""

	Lengths = proj['TrajLengths']

	sumInd = 0
	if len( data.shape ) > 1:
		if data.shape[1] ==1 :
			data = data[:,0]
		elif data.shape[0] ==1:
			data = data[0,:]

	data2d = -1 * ones( ( proj['NumTrajs'], Lengths.max() ) )

	for i in range( proj['NumTrajs'] ):
		data2d[ i, 0 : Lengths[i] ] = data[ sumInd : sumInd + Lengths[i] ]
		sumInd += Lengths[i]
	return data2d

def getDihedralFromAtoms( XYZ ):
	"""
	Calculate the dihedral angle between four atoms.

	Input:
	1) XYZ: np.array of shape ( 4, 3 ) containing four atoms' XYZ coordinates.
	
	Output:
	1) angle: angle (in radians) formed by the dihedral 1-2-3-4
	"""
	
	v12 = XYZ[0] - XYZ[1] # Vector from 2 -> 1
	v23 = XYZ[2] - XYZ[1] # Vector from 2 -> 3 ( could be from 3->2 this doesn't matter, we just use it to find perpindicular vectors )
	v34 = XYZ[3] - XYZ[2] # Vector from 3 -> 4
	
	# First project v12:
	#p12 = np.cross( v12, v23 ) # p12 is perpindicular to v12 and v23 
	p12 = [ v12[1]*v23[2] - v12[2]*v23[1], v12[2]*v23[0] - v12[0]*v23[2], v12[0]*v23[1] - v12[1]*v23[0] ]
	# Next, v34:
	#p34 = np.cross( v34, v23 ) # p34 is perpindicular to v34 and v23
	p34 = [ v34[1]*v23[2] - v34[2]*v23[1], v34[2]*v23[0] - v34[0]*v23[2], v34[0]*v23[1] - v34[1]*v23[0] ]
	# ^^^ Writing out explicitly, the cross products savest <2x time...
	#angle = np.arccos( np.dot( p34, p12 ) / np.sqrt( np.dot( p34, p34 ) * np.dot( p12, p12 ) ) ) # TIME: 3.130436e-04
	inside =  (p34[0]*p12[0] + p34[1]*p12[1] + p34[2]*p12[2] ) / np.sqrt( ( p34[0]*p34[0] + p34[1]*p34[1] + p34[2]*p34[2] ) * ( p12[0]*p12[0] + p12[1]*p12[1] + p12[2]*p12[2] ) )
	if ( abs( inside ) > 1 ):
		if abs( inside - ( - 1 ) ) < 1E-4: # i.e. inside ~ -1
			inside = -1
		elif abs( inside - ( 1 ) ) < 1E-4: # i.e. inside ~ 1
			inside = 1
		else:
			print "arccos( X ) is undefined because X = %e" % inside
	
	angle = np.arccos( inside ) 
	#angle = np.arccos( (p34[0]*p12[0] + p34[1]*p12[1] + p34[2]*p12[2] ) / np.sqrt( ( p34[0]*p34[0] + p34[1]*p34[1] + p34[2]*p34[2] ) * ( p12[0]*p12[0] + p12[1]*p12[1] + p12[2]*p12[2] ) ) ) # TIME : 2.050400e-04
	# ^^^ Writing out explicitly, the dot products saves 1.5x time...
	return angle
def AnalyzeFrameDihedrals( args ):
	frame, PsiAtomList, PhiAtomList, OmegaAtomList = args
	PsisFrame = [ getDihedralFromAtoms( frame[ ind ] ) for ind in PsiAtomList ]
	PsisFrame.insert( 0, -1 ) # Phis are not defined on the first amino acid
	PhisFrame = [ getDihedralFromAtoms( frame[ ind ] ) for ind in PhiAtomList ] 
	PhisFrame.append( -1 ) # Psis are not defined on the last amino acid
	OmegasFrame = [ getDihedralFromAtoms( frame[ ind ] ) for ind in OmegaAtomList ] 
	OmegasFrame.append( -1 ) # Omegas are not defined on the last amino acid
	#print PsisFrame, PhisFrame, OmegasFrame
	return PsisFrame, PhisFrame, OmegasFrame

def calcDihedrals( traj, procs = 1 ):
	""" 
	This function calculates the dihedral for a trajectory over time. The angle parameter gives which dihedral to calculate.

	Input:
	1) traj - A Trajectory.Trajectory object
	2) angle - One of [omega, psi, phi] indicating which angle to calculate

	Output:
	3) dihedrals - a np.array of shape len(traj), Nresidues
	"""

	atmNames = traj['AtomNames']
	backbone = np.where( ( atmNames == 'CA' ) + ( atmNames == 'C' ) + ( atmNames == 'N' ) )[0]

	traj.RestrictAtomIndices( backbone )

	nAtoms = traj['XYZList'].shape[1]
	# Now calculate all of the dihedrals
		# Phi : c - N - CA - C
		# Psi : N - CA - C - n
		# Omega : CA - C - n - ca
			# ^^^ NOTE: Capital letters mean on that particular amino acid, so notice that the Phi angle is 
			# NOT defined for the N-terminus, while the Psi and Omega angles are not defined for the C- terminus

	# Psi starts at N:
	PsiAtomList = [ np.arange( i, i+4 ) for i in np.where( traj['AtomNames'] == 'N' )[0] if i+4 <= nAtoms ]
	# Phi starts at C:
	PhiAtomList = [ np.arange( i, i+4 ) for i in np.where( traj['AtomNames'] == 'C' )[0] if i+4 <= nAtoms ]
	# Omega starts at CA:
	OmegaAtomList = [ np.arange( i, i+4 ) for i in np.where( traj['AtomNames'] == 'CA' )[0] if i+4 <= nAtoms ]
	# The above are the lists of indices to calculate the angles for.
	#print OmegaAtomList
	Psis = []
	Phis = []
	Omegas = []
	#print procs
	pool = mp.Pool( procs )
	N = traj['XYZList'].shape[0]
	result = pool.map_async( AnalyzeFrameDihedrals, zip( list( traj['XYZList'] ), N * [ PsiAtomList ], N * [ PhiAtomList ], N * [ OmegaAtomList ]  ) )
	pool.close()
	sol = result.get()
	pool.terminate()
	for item in sol:
		Psis.append( item[0] )
		Phis.append( item[1] )
		Omegas.append( item[2] )	

	
	return array( Psis ), array( Phis ), array( Omegas )

def getContactMapHelper( args ):
	"""
	This is a helper function for using multiprocessing and getContactMap
	"""
	XYZ = args[0]
	resList = args[1]
	cutoff = args[2] 

	return getContactMap( XYZ, resList = resList, cutoff = cutoff )

def getContactMap( XYZ, resList = None, cutoff = 0.6 ):
	"""
	This function calculates the contact map for a given XYZ coordinate of atoms. NOTE: This calculates contacts between atoms unless resList is given which specifies to calculate residue contacts by the closest distance between atoms of any two residues.

	Input:
	1) XYZ - np.array of xyz coordinates for the atoms in question ( units = nm )
	2) resList [ None ] - list of np.arrays of indices in XYZ corresponding to which atoms correspond to which residue. So the output matrix will be N x N where N = len( resList ).
	3) cutoff [ 0.6 ] - Cutoff to use in calculating the contacts ( units = nm )

	Output:
	3) ContactMap - np.array of size N x N, which contains only the upper triangular values in the contact map minus the diagonal, and 2 off-diagonals (corresponding to residues 3 places away from eachother)
	"""
	cutoff2 = cutoff * cutoff # Avoid the sqrt, and compare to cutoff**2
	# If no resList, then assume each atom represents a residue
	if resList == None:
		resList = np.array( [ [i] for i in range( len( XYZ ) ) ] )
  
	PW_DistAtoms = np.ones( ( len( XYZ ), len( XYZ ) ) ) * cutoff2 
	Inter = np.zeros( ( len( resList ), len( XYZ ) ) ) 
	PW_Dist = np.zeros( ( len( resList ), len( resList ) ) )  
	# In the end will set to zero anything that is BIGGER than the cutoff, so the diagonals should be set to the cutoff2 so that they are always set to 1
   
	# We will determine the contacts by computing the pairwise distances:
	for ind, atom in enumerate( XYZ[:-len( resList[-1] ) ] ): # Go up to the last residue's atoms
	# First find out which residue this atom is in:
		inRes = np.where( np.array( [ ( ind in res ) for res in resList ] ) )[0]
		if len( inRes ) != 1:
			print "atom %d in multiple residues... Exiting..."
			exit()
		inRes = inRes[0]

		# Since we only need the upper triangular, three diagonals we need only calculate the distance to a portion of the other atoms.
		#print resList[ inRes + 1 : ]
		#print XYZ[ resList[ inRes + 1 : ] ]
		TestInd = np.concatenate( resList[ inRes + 1 : ] ) # These are the indices of the residues greater than this index's residue
		TestAtoms = np.row_stack( XYZ[ TestInd ] ) # These are the coordinates for the aboe residues
		diff = TestAtoms - atom
		dists = ( diff * diff ).sum( axis = 1 )
		
		PW_DistAtoms[ ind, TestInd ] = dists
		PW_DistAtoms[ TestInd, ind ] = dists
		
	PW_DistAtoms = (PW_DistAtoms <= cutoff2)

	# since the atoms belong to residues, we have to turn each block into a signle value based on whether that block as all zeros or at least one one. If it has a one in it, then the residues are contacting
	for i in range( len( resList ) ):
		Inter[ i, : ] = PW_DistAtoms[ resList[i] , : ].any( axis = 0 ) # First contract axis = 0 ( rows )
	for i in range( len( resList ) ):
		PW_Dist[ i, : ] = Inter[ :, resList[i] ].any( axis = 1 ) # Then contract axis = 1 ( columns )
	return PW_Dist

def calcContacts( traj, procs = 1, cutoff = 0.6, RestrictToCA=False):
	"""
	This function calculates the contact map for an entire trajectory, and returns it as a np.array.

	Inputs:
	1) traj - Trajectory.Trajectory object
	2) procs [ 1 ] - number of processes to run using python.multiprocessing
	3) cutoff [ 0.6 ] - cutoff to define a contact, NOTE: this cutoff is the distance between any two atoms in a trajectory.
	4) RestrictToCA [ False ] - If true then the trajectory will be restricted to only CA's which means there is only one distance to calculate between residues. This will be faster, but since some residues are different sizes, a constant cutoff doesn't really make sense...

	Outputs:
	1) contactMaps - np.array of dimension ( N, m, m ) where N = number of frames in traj, and m is the number of residues in traj 
	"""

	pool = mp.Pool( procs )

	if RestrictToCA:
		traj.RestrictAtomIndices( np.where( traj['AtomNames'] == 'CA' )[0] )

	# Need to produce a residue list corresponding to the atom indices in each residue for this to work.
	# Doesn't matter how the residues are numbered, as long as they are IN ORDER... this should be a safe assumption...
	resIDs = traj['ResidueID']
	if ( ( resIDs[1:] - resIDs[:-1] ) < 0 ).any():
		print "Residue IDs in trajectory are NOT in order! Fix this and then call this function again. You shouldn't have miss-ordered residue indices anyway..."
		exit()

	uniqueResIDs = np.unique( resIDs ) # This SORTs the id's so if they are not in order, then we are in trouble...
	resList = [ np.where( resIDs == i )[0] for i in uniqueResIDs ] # Get the indices for each residue

	print traj['XYZList'].shape
	N = traj['XYZList'].shape[0]

	#result = pool.map_async( getContactMap, traj['XYZList'] )
	result = pool.map_async( getContactMapHelper, zip( traj['XYZList'], N * [ resList ], N * [ cutoff ] ) )
	pool.close()
	pool.join()
	sol = result.get()

	contactMaps = np.array( sol )

	return contactMaps

def getAtomDistanceFromFrame( frame ):
	""" 
	This function returns the atomic distances from a single frame
	"""
	n1 = frame.shape[0]
	temp = np.zeros( ( n1, n1, 3 ) )
	
	for row in xrange( n1 ):
		temp[ row, row+1 :, : ] = frame[ row ] - frame[ row+1 : ]
	temp *= temp
	temp = temp.sum(axis=2)

	return temp

def getAllAtomDistances( XYZList, nProcs = 1 ):
	"""
	This function returns all intramolecular atomic distances (L2 norm in 3D) as a matrix. To save space, only the upper triangular component is calculated.

	Inputs:
	1) XYZList - np.array with a trajectory's list of frames. The shape is (n0, n1, 3), where n0 is the number of frames, and n1 is the number of atoms.

	Outputs:
	2) AA_Traj - np.array with a trajectory's all atom distances calculated. The shape is (n0, n1, n1), where n0 and n1 are as above. NOTE: The matrix is upper triangular!
	"""

	n0,n1,n2 = XYZList.shape

	pool = mp.Pool( nProcs )
	result = pool.map_async( getAtomDistanceFromFrame, XYZList )
	sol =result.get()
	pool.close()
	AA_Traj = np.concatenate( sol )

	a=time()
	AA_Traj = np.sqrt( AA_Traj )
	print time()-a
	del sol
	return AA_Traj

def distributeStateAvgs( data, pops, bins=200, xLimit=None ):
   """
   This function takes in an Nx3 array which has the state number, average value, standard deviation and adds gaussians corresponding to these values with the indicated integral in pops

   Inputs:
   1) data - Nx3 array of data values should be state number (not used), average value, standard devaition
   2) pops - Nx1 array of populations for each state to use for what each gaussian integrates to.

   Output:
   1) distXY - 200 x 2 array of data points corresponding to x values and y values which are the distribution heights
   """

   means = data[:,1] 
   stds = data[:,2]
   
   # Limit to where stds != 0:
   goodInd = np.where( stds != 0 )
   means = means[ goodInd ]
   stds = stds[ goodInd ]
   pops = pops[ goodInd ]

   if len( pops.shape ) != 1:
      pops = pops[:,0]

   pops = pops / ( np.sqrt( 2. * np.pi ) * stds )
   if not xLimit:
      minX = np.max( ( 0, ( means - 2 * stds ).min() ) )
      maxX = np.max( means + 2 * stds )
   else:
      minX, maxX = xLimit

   x = np.linspace( minX, maxX, bins )
   y = np.zeros( len( x ) )

   for i in xrange( len( means ) ):
      y +=  FitData.gauss( [ pops[i], means[i], stds[i] ], x )

   return np.array( zip(x,y) )

def calc_time_avg(traj, chunk_size ):
   traj_shape = traj.shape
   N=traj_shape[0]
   if not N % chunk_size:
      num_chunks = N / chunk_size
   else:
      num_chunks = int( N / chunk_size + 1 )
   avg_traj_shape = list( traj_shape )
   avg_traj_shape[0] = num_chunks
   avg_traj = np.zeros( avg_traj_shape ).astype( traj.dtype )
   for i in range(num_chunks):
       chunk = traj[ i * chunk_size : ( i + 1 ) * chunk_size ]
       chunk = chunk.mean(axis=0)
       avg_traj[i] = chunk

   return avg_traj

