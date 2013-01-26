#!/usr/bin/env python 

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-P',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option('-a',dest='ass_FN',default='Assignments.Fixed.h5',help='Assignments from msmbuilder [ ./Assignments.Fixed.h5 ]')
parser.add_option('-x',dest='x_FN',help='X data for diffusion (Raw data)')
parser.add_option('-y',dest='y_FN',help='Y data for diffusion (Raw data)')
parser.add_option('-p',dest='top_FN',help='Topology to use with running simulations')
parser.add_option('-f',dest='mdp_FN',help='mdp file to use with running simulations')
parser.add_option('-t',dest='table_FN',default='./table.xvg',help='Table for running CAlpha simulations. Pass None if you don\'t want to use a table. [ ./table.xvg ]')
parser.add_option('--co',dest='cutoff_FN',help='Cutoffs for each native contact')
parser.add_option('--nc',dest='nc_FN',help='Native Contacts file from Onuchic webserver')
parser.add_option('--rr',dest='ranges',action='append',help='Residue ranges for defining folded/unfolded')
parser.add_option('--fc',dest='fc_FN',default='./FCs.dat',help='Forward committors from msmbuilder')
parser.add_option('--fx',dest='fold_x',type=float,help='Cutoff for folded state in x-axis')
parser.add_option('--ux',dest='unfold_x',type=float,help='Cutoff for unfolded state in x-axis')
parser.add_option('--fy',dest='fold_y',type=float,help='Cutoff for folded state in y-axis')
parser.add_option('--uy',dest='unfold_y',type=float,help='Cutoff for unfolded state in y-axis')
parser.add_option('--low-is-folded',dest='low_is_folded',default=False,action='store_true',help='Pass this flag if your reaction coordinate is low when folded (eg. RMSD to Native state)')
parser.add_option('-s',dest='state',type=int,help='State in msm to calculate pfolds for')
options, args = parser.parse_args()

import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
import os, sys, re
from msmbuilder import Serializer, Project
import subprocess
from pyschwancr import dataIO, MonteCarlo

os.putenv('GMX_MAXBACKUP','-1')

def SetupMover( xDat, yDat ):
	# Set up the mover:
	MC_Mover = MonteCarlo.MC_Mover( xDat, yDat )
	# Determine folded and unfolded state indices.
	NoX = False
	NoY = False
	if options.fold_x and options.unfold_x: # If  we defined both on the command line then use them
		if options.low_is_folded: # If low is folded then we want:
			X_f = np.array( np.where( MC_Mover.xTicks <= options.fold_x )[0] ) # All indices below the folded cutoff
			X_u = np.array( np.where( MC_Mover.xTicks >= options.unfold_x )[0] ) # All indices above the unfolded cutoff
		else: # If low is not folded
			X_f = np.array( np.where( MC_Mover.xTicks >= options.fold_x )[0] ) # All the indices above the folded cutoff
			X_u = np.array( np.where( MC_Mover.xTicks <= options.unfold_x )[0] ) # All the indices below  the unfolded cutoff
	else:
		NoX = True # If we didn't use the command line arguments, then pick all of the indices
		X_f = np.arange( len( MC_Mover.xTicks ) )
		X_u = np.arange( len( MC_Mover.xTicks ) )
	# Repeat for the y-axis
	if options.fold_y and options.unfold_y:
		if options.low_is_folded:
			Y_f = np.array( np.where( MC_Mover.yTicks <= options.fold_y )[0] )
			Y_u = np.array( np.where( MC_Mover.yTicks >= options.unfold_y )[0] )
		else:
			Y_f = np.array( np.where( MC_Mover.yTicks >= options.fold_y )[0] )
			Y_u = np.array( np.where( MC_Mover.yTicks <= options.unfold_y )[0] )
	else:
		NoY = True
		Y_f = np.arange( len( MC_Mover.yTicks ) )
		Y_u = np.arange( len( MC_Mover.yTicks ) )

	if NoX and NoY: # If neither of the x or y axes were defined for 
		print "Need to define folded and unfolded states along at least one axis!"
		exit()
	return MC_Mover, X_f, X_u, Y_f, Y_u

def MC_Pfold( MC_Mover, r0, X_f, X_u, Y_f, Y_u, plotLast = False ):

	Pfold = MonteCarlo.CalcPfold( MC_Mover, r0, X_f, X_u, Y_f, Y_u, Nsteps = 2000, Ntrajs = 100 )

	if plotLast:
		imshow( np.log10( MC_Mover.Pot.T ), origin='lower', extent = [0,1,0,1])
		traj = MC_Mover.getXYvalues()
		plot( traj[:,0], traj[:,1], color='white',marker='o' )
		vlines( options.fold_x, 0, 1, color='black')
		vlines( options.unfold_x, 0, 1, color='black')
		savefig('LastTraj.pdf')

	return Pfold

def StartRawTrajs( whichConfs, Proj ):
	# First make the dnirectory to save things.
	if not os.path.exists( './RawTrajs' ):
		os.system( 'mkdir ./RawTrajs' )
	Pfolds = []
	count = 0.
	for trajFrame in whichConfs:
		print "traj: %d, frame: %d" % tuple( trajFrame ),
		tempTraj = Proj.GetConformations( trajFrame.reshape( (1,2 ) ) )
		if not os.path.exists(  os.path.join( 'RawTrajs', 'State%d_trj%d_frm%d' % ( options.state, trajFrame[0], trajFrame[1] ) ) ):
			os.system( "mkdir %s" % ( os.path.join( 'RawTrajs', 'State%d_trj%d_frm%d' % ( options.state, trajFrame[0], trajFrame[1] ) ) ) )
		tempTraj.SaveToPDB( 'RawTrajs/State%d_trj%d_frm%d/conf.pdb' % ( options.state, trajFrame[0], trajFrame[1] ) )
		cmd = "name=State%d_trj%d_frm%d;" % ( options.state, trajFrame[0], trajFrame[1] ) + "editconf -f RawTrajs/$name/conf.pdb -o RawTrajs/$name/conf.gro -c -box 11.31800  11.22900  11.30200 &> /dev/null" ##### Need to edit this box for different systems. This is for CheY.
		os.system( cmd )
	#	print "Made PDBs...",
		# ^^^^^^^ Convert the pdb's to gro files with the correct box
		cmd = "name=State%d_trj%d_frm%d;" % ( options.state, trajFrame[0], trajFrame[1] ) + " for i in {0..99}; do grompp -f %s -c ./RawTrajs/$name/conf.gro -p %s -maxwarn 1 -o ./RawTrajs/$name/$i.tpr &> /dev/null; done"  % ( options.mdp_FN, options.top_FN )
		os.system( cmd )
		# ^^^^^^ Make the tpr files.
	#	print "Converted to GROs",
		#print "Converted PDBs to GROs ... with box: %s" % ("11.31800  11.22900  11.30200 " ),
		# Now submit these jobs to run in a pbs queue:
		filename = 'State%d_trj%d_frm%d' % ( options.state, trajFrame[0], trajFrame[1] )
		dirname =os.path.join( os.path.abspath( '.' ), 'RawTrajs', filename )
		cmd = "cp %s %s" % ( options.table_FN, dirname )
		os.system( cmd )
		
		print "Running simulations", 
		cmd = "runSims.sh %s 0 100 &> /dev/null" % ( dirname ) 
		os.system( cmd )

	#	print "Calculating Raw pfold.",
		pfoldOpts = '-t %s --uc %s --fc %s --nc %s --co %s -P 24 --rr %s --rr %s' % ( dirname, options.unfold_x, options.fold_x, options.nc_FN, options.cutoff_FN, options.ranges[0], options.ranges[1] )
	
		cmd = "calcPfold_AnyData.py " + pfoldOpts

		pfoldOutput = subprocess.check_output( cmd, shell=True )
	
		matchObj = re.search( 'Folded\s=\s(\d+)\sUnfolded\s=\s(\d+)', pfoldOutput )
		u = int( matchObj.group(1) )
		f = int( matchObj.group(2) )
		if u+f < 85:
			print "Less than 85% of the raw trajectories ended. You should increase nsteps"
		Pfold = float(u) / float( u + f )
		Pfolds.append( Pfold )
		count += 1. / len( whichConfs )
		print "Calculating Raw Pfolds. %.2f%% Done.\r" % ( count * 100 ),

	return Pfolds
	
def main():

	# First load all the data.
	print "Loading data ...",
	Proj = Project.Project.LoadFromHDF( options.proj_FN )
	Ass = Serializer.LoadData( options.ass_FN )
	xDat = dataIO.readData( options.x_FN )[:,0]
	yDat = dataIO.readData( options.y_FN )[:,0]
	print "Done."

	print "Reformatting raw data...",
	xFmtd = np.ones( Ass.shape ) * -1
	yFmtd = np.ones( Ass.shape ) * -1
		
	lengthSum = 0
	for i in range( len( Proj['TrajLengths'] ) ):
		#print xDat[ lengthSum : lengthSum + Proj['TrajLengths'][i] ].shape
		#print xFmtd[i][:Proj['TrajLengths'][i]].shape
		xFmtd[i,:Proj['TrajLengths'][i]] = xDat[ lengthSum : lengthSum + Proj['TrajLengths'][i] ]
		yFmtd[i,:Proj['TrajLengths'][i]] = yDat[ lengthSum : lengthSum + Proj['TrajLengths'][i] ]
		lengthSum += Proj['TrajLengths'][i]
	print "Done."

	MC_Mover, X_f, X_u, Y_f, Y_u = SetupMover(xDat, yDat)
	if options.state:
		whichState = options.state
	else:
		whichState = np.random.randint( Ass.max() + 1 )

	whichConfs = np.array( np.where( Ass == whichState ) ).T
	if whichConfs.shape[0] > 100:
		whichConfs = whichConfs[:100]
		print "Only using the first 100 conformations in the state... (Out of %d)" % ( len( np.where( Ass == whichState )[0] ) )
	# whichConfs will be used to get the pdb's for each conformation.
	Pfolds_Raw = StartRawTrajs( whichConfs, Proj )

	startingPositionsXY = np.array( zip( xFmtd[ np.where( Ass == whichState )	], yFmtd[ np.where( Ass == whichState ) ] ) )
	# These are the starting locations for doing diffusion.
	Pfolds_MC = []
	count = 0.
	for XY in startingPositionsXY:
		count += 1. / len( startingPositionsXY )
		print "Calculating Diffusion Pfolds. %.2f%% Done.\r" % ( count * 100 ),
		r0 = MC_Mover.getIndices( XY[0], XY[1] )
		Pfolds_MC.append( MC_Pfold( MC_Mover, r0, X_f, X_u, Y_f, Y_u ) )

	xDat_name = '.'.join( options.x_FN.split('/')[-1].split('.')[:-1] )
	yDat_name = '.'.join( options.y_FN.split('/')[-1].split('.')[:-1] )
	np.savetxt('MC_Pfolds_State%d_%s_%s.dat' % ( whichState, xDat_name, yDat_name ), Pfolds_MC )
	print "Saved MC_Pfolds to %s" % 'MC_Pfolds_State%d_%s_%s.dat' % ( whichState, xDat_name, yDat_name )
	np.savetxt('Raw_Pfolds_State%d.dat' % options.state, Pfolds_Raw )
	print "Saved Raw_Pfolds to %s" % 'Raw_Pfolds_State%d.dat' % whichState 

	plot( Pfolds_MC, Pfolds_Raw, '.' )
	xlabel('Diffusion Pfolds')
	ylabel('Raw Pfolds (MD)')
	title('Raw Pfolds vs Diffusion Pfolds')
	xlim([0,1])
	ylim([0,1])
	savefig('State%d_RawVsDiff_%s_%s.pdf' % ( whichState, xDat_name, yDat_name ) )
	print "Saved plot to %s" % 'State%d_RawVsDiff_%s_%s.pdf' % ( whichState, xDat_name, yDat_name )



if __name__ == '__main__':
	main()


