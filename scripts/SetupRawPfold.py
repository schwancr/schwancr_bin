#!/home/schwancr/epd-7.1-1-rh5-x86_64/bin/python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option('-a',dest='ass_FN',default='./Assignments.Fixed.h5',help='Assignments from msmbuilder [ ./Assignments.Fixed.h5 ]')
parser.add_option('-s',dest='states',type=int,action='append',help='State to sample from the msm. You may pass multiple -s\'s')
parser.add_option('-c',dest='num_confs',type=int,default=100,help='Number of conformations to sample [ 100 ]')
parser.add_option('-b',dest='gro_FN',help='Gro file to use in determining the box size. The box size is the only thing that is used, so the protein in it doesn\'t matter')
parser.add_option('--top',dest='top',help='Top file for grompp')
parser.add_option('-f',dest='mdp',help='mdp file for grompp')
parser.add_option('-t',dest='table',help='table.xvg if your sims need it. WILL CREATE SYMLINKS TO THIS FILE')
options, args = parser.parse_args()

from msmbuilder import Serializer, Project
from numpy import *
import os, sys, re
 
os.putenv( 'GMX_MAXBACKUP','-1' )
Ass = Serializer.LoadData( options.ass_FN )
Proj = Project.Project.LoadFromHDF( options.proj_FN )

boxString = open( options.gro_FN, 'r' ).readlines()[-1].strip()

for state in options.states:
	print "Working on state %6d\r" % state,
	try: os.mkdir( 'State%d' % state )
	except: 
		print "State%d already exists. Continuing" % state
		continue
	os.mkdir( 'State%d/Confs' % state )
	os.mkdir( 'State%d/Sims' % state )
	os.mkdir( 'State%d/PDBs' % state )

	which = array( where( Ass == state ) ).T
	if len( which ) < options.num_confs:
		print "Only %d conformations in state %d" % ( options.num_confs, state )
	else:
		which = random.permutation( which )[:options.num_confs] # Random sample of the conformations
	
	count = 0
	for pair in which:
		count += 1
		print 'Working on state %6d, conformation %3d\r' % (state, count),
		conf = Proj.GetConformations( array([ pair ]) )
		pdbName = 'State%d/PDBs/State%d_trj%d_frm%d.pdb' % ( state, state, pair[0], pair[1] )
		conf.SaveToPDB( pdbName )
		groName = 'State%d/Confs/State%d_trj%d_frm%d.gro' % ( state, state, pair[0], pair[1] )
		
		os.system('editconf -c -f %s -o %s -box %s &> /dev/null'  % ( pdbName, groName, boxString ) )

		dirName = 'State%d/Sims/trj%d_frm%d' % (state, pair[0], pair[1] )
		os.mkdir( dirName )
	
		os.system('for i in {0..99}; do grompp -f %s -c %s -p %s -o %s/$i.tpr -maxwarn 1 &> /dev/null; done' % ( options.mdp, groName, options.top, dirName ) )

		if options.table:
			os.system('ln -s %s %s' % ( options.table, dirName ) )
	print "\n"	
	print "Submitting directories ..."
	os.system('submitPfolds.sh State%d' % (state) )
