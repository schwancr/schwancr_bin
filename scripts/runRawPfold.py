#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-P',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option('-n',dest='num_confs',default=100,type=int,help='Number of conformations to sample from the project [ 100 ]')
parser.add_option('-f',dest='mdp_FN',default='./pfold.mdp',help='MDP to run the simulations with [ ./pfold.mdp ]')
parser.add_option('-p',dest='top_FN',help='Topology file for gromacs')
parser.add_option('-t',dest='table_FN',default='./table.xvg',help='Table for running sims. [ ./table.xvg ]')
parser.add_option('-c',dest='gro_FN',default='./conf.gro',help='Sample gro file. ONLY USED FOR BOX SIZE [ ./conf.gro ]' )
options, args = parser.parse_args()
 
from numpy import *
from msmbuilder import Project
from pyschwancr import dataIO
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

if not os.path.exists( 'RawTrajs' ):
	os.mkdir( 'RawTrajs' )
totalConfs = 0
while ( totalConfs < options.num_confs ):
# Loop until we've actuall done 100 conformations
	# First grab a random trajectory, then a frame
	
	trajNum = random.randint( len( Proj['TrajLengths'] ) )
	frameNum = random.randint( Proj['TrajLengths'][ trajNum ] )
 
	print "Trying trajectory %d frame %d" % ( trajNum, frameNum )
	dirName = 'trj%d_frm%d' % ( trajNum, frameNum )

	for tempDir in dirs2check:
		testAry = array( [ re.search( dirName, thing ) for thing in os.listdir( tempDir ) ] )
		if testAry.any():
			cmd = 'ln -s %s/*_%s `pwd`/RawTrajs/%s' % ( tempDir, dirName, dirName )
			os.system( cmd )
			totalConfs += 1
			print "Found directory! This conf was already simulated. Linking the directory and going to the next conf"
			break

	# Check to see if this conf has already been calculated in this directory:
	if os.path.exists( os.path.join( 'RawTrajs',dirName) ):
		print "Trajectory %d and Frame %d has already been done. Continuing..."
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
	cmd = "for i in {0..99}; do grompp -f %s -c RawTrajs/%s/conf.gro -p %s -maxwarn 1 -o RawTrajs/%s/$i.tpr &> /dev/null & done" % ( options.mdp_FN, dirName, options.top_FN, dirName )
	os.system( cmd )

	# Make a link to the table file so runSims.sh can find it.
	cmd = "ln -s `pwd`/%s `pwd`/RawTrajs/%s/" % ( options.table_FN, dirName )
	os.system( cmd )

	# Run the simulations:
	cmd = "runSims.sh $(pwd)/RawTrajs/%s/ 0 100" % dirName
	os.system( cmd )	

	totalConfs += 1



