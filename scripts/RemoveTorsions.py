#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p',dest='top_FN',help='Topology file name to edit')
parser.add_option('-c',dest='gro_FN',help='Gro file name to use to calculate secondary structure')
parser.add_option('-s',dest='SS_list',action='append',help='List of secondary structure to KEEP torsions/angles for')
parser.add_option('-o',dest='out_FN',help='Output file name')
options, args = parser.parse_args()
 
import os, sys, re
 
top = open( options.top_FN, 'r' )

if options.out_FN:
	outTop = open( options.out_FN, 'w' )
else:
	outTop = open( options.top_FN[:-4] + '_keep_' + '-'.join( options.SS_list ) + '.top' , 'w' )
# First generate the ss assignments:
os.putenv('GMX_MAXBACKUP','-1')
cmd = "echo \"Protein\" | do_dssp -f %s -s %s " % ( options.gro_FN, options.gro_FN )

os.system( cmd )

ss = open( 'ss.xpm', 'r' )

SS_str = ''
 #Read in the list of assignments:
for line in ss:
	matchObj = re.search( '^\"([SHT~])\"', line )

	if matchObj: # Then we have a residue!!!
		SS_str += matchObj.group(1)
		
ss.close()

ResToKeep = []
for i,ass in enumerate(SS_str):
	if ass in options.SS_list:
		#print ass, i + 1
		ResToKeep.append( i + 1 )
	#else:
		#print ass, i + 1, '*'


for line in top:
	outTop.write(line)
	if re.search( '\[ angles \]', line ):
		break
	
# In angles:
for line in top:
	if re.search( '\[ dihedrals \]', line ):
		outTop.write(line)
		break
	elif line.strip() == '':
		outTop.write(line)
		continue
	else:
		matchObj = re.search( '^\s*(\d+)\s*(\d+)\s*(\d+)', line )
		if matchObj:
			atoms = [ int( i ) for i in matchObj.groups() ]
		
			keepAtoms = [ i in ResToKeep for i in atoms ]

			if sum( keepAtoms ) == 3:
				outTop.write( line )
			
# In dihedrals
for line in top:
	matchObj = re.search( '^\s*(\d+)\s*(\d+)\s*(\d+)\s*(\d+)', line )
	if matchObj:
		atoms = [ int( i ) for i in matchObj.groups() ]
		
		keepAtoms = [ i in ResToKeep for i in atoms ]

		if sum( keepAtoms ) == 4:
			outTop.write( line )
	else:
		outTop.write( line )

outTop.close()
top.close()
