#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-P',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option('-f',dest='mdp_FN',default='./pfold.mdp',help='MDP to run the simulations with [ ./pfold.mdp ]')
parser.add_option('-p',dest='top_FN',help='Topology file for gromacs')
parser.add_option('-t',dest='table_FN',default='./table.xvg',help='Table for running sims. [ ./table.xvg ]')
parser.add_option('-c',dest='gro_FN',default='./conf.gro',help='Sample gro file. ONLY USED FOR BOX SIZE [ ./conf.gro ]' )
parser.add_option('-x',dest='x_dat',help='X-data to use in filling the grid' )
parser.add_option('-y',dest='y_dat',help='Y-data to use in filling the grid' )
parser.add_option('--x-interval',dest='x_int',help='Interval to look at x-data ( --x-interval a,b )' )
parser.add_option('--y-interval',dest='y_int',help='Interval to look at y-data ( --y-interval v,w )' )
options, args = parser.parse_args()
 
from numpy import *
from msmbuilder import Project
from pyschwancr import dataIO, msmTools
import os, sys, re
 

Proj = Project.Project.LoadFromHDF( options.proj_FN )

groFile = open( options.gro_FN, 'r' )
boxLine = groFile.readlines()[-1].strip()
print "Using box %s" % boxLine
groFile.close()

os.putenv( 'GMX_MAXBACKUP', '-1' )

dirs2check = []
dirs2check.append( '/home/schwancr/CheY_smog/fahCheY/Calpha/SkipFrame1/Data_1.5A_u100/RawTrajs' )
dirs2check.append( '/home/schwancr/CheY_smog/fahCheY/Calpha/SkipFrame1/Data_1.5A_u100/5_5_TPT/State24814_SIMS/Sims' )
dirs2check.append( '/home/schwancr/CheY_smog/fahCheY/Calpha/SkipFrame1/Data_1.5A_u100/5_5_TPT/State25033/Sims' )
dirs2check.append( '/home/schwancr/CheY_smog/fahCheY/Calpha/SkipFrame1/RawVsDiffusion/Data_8A_u10/regAssign/RawTrajs' )
dirs2check.append( '/home/schwancr/CheY_smog/fahCheY/Calpha/SkipFrame1/RawVsDiffusion/Data_8A_u10/regAssign/State92830/RawTrajs' )
dirs2check.append( '/home/schwancr/CheY_smog/fahCheY/Calpha/SkipFrame1/RawVsDiffusion/RandomSampling/RawTrajs' )

if not os.path.exists( 'RawTrajs' ):
	os.mkdir( 'RawTrajs' )


def run( TrajFrame ):
	totalConfs = 0
	for pair in TrajFrame:
	
		trajNum = pair[0]
		frameNum = pair[1] 
 
		print "Trying trajectory %d frame %d" % ( trajNum, frameNum )
		dirName = 'trj%d_frm%d' % ( trajNum, frameNum )

		for tempDir in dirs2check:
			testAry = array( [ re.search( dirName, thing ) for thing in os.listdir( tempDir ) ] )
			if testAry.any():
				cmd = 'ln -s %s/*%s `pwd`/RawTrajs/%s' % ( tempDir, dirName, dirName )
				os.system( cmd )
				totalConfs += 1
				print "Found directory! This conf was already simulated. Linking the directory and going to the next conf"
				break

	# Check to see if this conf has already been calculated in this directory:
		if os.path.exists( os.path.join( 'RawTrajs',dirName) ):
			print "Trajectory %d and Frame %d has already been done. Continuing..." % ( pair[0], pair[1] )
			continue

		os.mkdir( os.path.join( 'RawTrajs', dirName ) )

	# Save the pdb:
		Traj = Proj.GetConformations( array([ [ trajNum, frameNum ] ]) )
		Traj.SaveToPDB( os.path.join( 'RawTrajs', dirName, 'conf.pdb' ) )
	
	# Add the box:
		cmd = "editconf -f RawTrajs/%s/conf.pdb " % dirName
		cmd += "-c -box %s " % boxLine
		cmd += "-o RawTrajs/%s/conf.gro &>/dev/null &" % dirName 
		os.system( cmd )

	# Now generate the tprs
		cmd = "for i in {0..99}; do grompp -f %s -c RawTrajs/%s/conf.gro -p %s -maxwarn 1 -o RawTrajs/%s/$i.tpr &> /dev/null; done" % ( options.mdp_FN, dirName, options.top_FN, dirName )
		os.system( cmd )

	# Make a link to the table file so runSims.sh can find it.
		cmd = "ln -s `pwd`/%s `pwd`/RawTrajs/%s/" % ( options.table_FN, dirName )
		os.system( cmd )

	# Run the simulations:
		cmd = "runSims.sh $(pwd)/RawTrajs/%s/ 0 100" % dirName
		os.system( cmd )	

		totalConfs += 1

def main():
	# Need to construct the list of trajframes to use. First look at what the values are for the ones already done.

	OldTrajFrames = []
	
	for tempDir in dirs2check:
		tempList = [ fn for fn in os.listdir( tempDir ) if re.search( 'trj\d+_frm\d+$', fn ) ]
		tempTrajFrames = [ re.search( 'trj(\d+)_frm(\d+)$', fn ).groups() for fn in tempList ]
		tempTrajFrames = [ [ int( a ), int( b ) ] for (a,b) in tempTrajFrames ]

		OldTrajFrames.extend( tempTrajFrames )

	Xdat = dataIO.readData( options.x_dat )
	Ydat = dataIO.readData( options.y_dat )
	uniqX = unique( Xdat )
	uniqY = unique( Ydat )
	diffX = abs( uniqX[1:] - uniqX[:-1] ).min()
	diffY = abs( uniqY[1:] - uniqY[:-1] ).min()
	
	Nx = 1. / diffX
	Ny = 1. / diffY

	Xdat = Xdat * ( Nx ) 
	Ydat = Ydat * ( Ny ) 
	
	if len( Xdat.shape ) > 1:
		Xdat = Xdat[:,0]
	if len( Ydat.shape ) > 1:
		Ydat = Ydat[:,0]

	Xdat2D = msmTools.reshapeRawData( Xdat.astype(int), Proj )
	Ydat2D = msmTools.reshapeRawData( Ydat.astype(int), Proj )
	
	x_interval = [ float( i ) for i in options.x_int.split(',') ]
	y_interval = [ float( i ) for i in options.y_int.split(',') ]

	x_range = arange( int( x_interval[0] * Nx ), int( x_interval[1] * Nx ) + 1 )
	y_range = arange( int( y_interval[0] * Ny ), int( y_interval[1] * Ny ) + 1 )

	#print x_range, y_range
	#print Xdat.max(), Ydat.max()
	#print OldTrajFrames
	# Now construct the list of traj frames to use in the analysis
	TrajFrames = []
	print "Finding conformations..."
	for Xi in x_range:
		for Yi in y_range:
			whichTrajFrames = array( where( ( Xdat2D == Xi ) * ( Ydat2D == Yi) ) ).T
			if len( whichTrajFrames ) <= 4:
				ToAddTrajFrames = [ list( i ) for i in whichTrajFrames ] 
			else:
				ToAddTrajFrames = [ list( pair ) for pair in whichTrajFrames if list( pair ) in OldTrajFrames ] 
				# The above list contains pairs for this x,y for which the simulations have already been done.
				if len( ToAddTrajFrames ) > 4:
					ToAddTrajFrames = [ list( pair ) for pair in random.permutation( ToAddTrajFrames )[:4] ]
				while len( ToAddTrajFrames ) < 4:
					randPair = whichTrajFrames[ random.randint( len( whichTrajFrames ) ) ]
					randPair = list( randPair )
		#			print '\t', randPair, ToAddTrajFrames
					if randPair in ToAddTrajFrames:
						continue
					else:
						ToAddTrajFrames.append( randPair )
			TrajFrames.extend( ToAddTrajFrames )				
	#		print Xi, Yi, ToAddTrajFrames, [ ( Xdat2D[ tuple(pair) ], Ydat2D[ tuple(pair) ] ) for pair in ToAddTrajFrames ]

	print TrajFrames

	print "Running the simulations..."
	run( TrajFrames )

if __name__ == '__main__':
	main()

