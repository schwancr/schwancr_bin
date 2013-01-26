#!/usr/bin/env python

import sys, os, re

inputFile = open( sys.argv[1], 'r' )

lines = inputFile.readlines()

# First generate the dictionary, which is a dictionary with residue dictionaries.

RenameDict = {}
CatchAllDict = { 'OXT': 'O2' }
inAtoms = False
tempDict = {}
currentRes = ''
for line in lines:
	if re.search( '^;',line ):
		continue
	elif re.search( '^\[ bondedtypes \]', line ):
		continue
	elif line.strip() == '':
		continue
	elif re.search( '\s+\d+\s+\d+', line ):
		continue

	resMatch = re.search( '^\[\s(.*)\s\]', line )
	if resMatch:
		if currentRes:
			RenameDict[ currentRes ] = tempDict
			tempDict = {}
			
		currentRes = resMatch.group(1).strip()
		inAtoms = False

	inAtmsMatch = re.search( '^\s\[\satoms\s\]', line )
	if inAtmsMatch:
		inAtoms = True
		continue

	inBondsMatch = re.search( '^\s\[\sbonds\s\]', line )
	if inBondsMatch:
		inAtoms = False

	if inAtoms:
		tempDict[ line.split()[0] ] = line.split()[1]
	else:
		continue
		
# Now the dictionary is made. Use it to rewrite the topology in sys.argv[2]

top = open( sys.argv[2], 'r' )
outTop = open( sys.argv[2][:-4] + '_NEW.top' , 'w' )

inAtomTypes = False
for line in top:
	
	if re.search( '^\s*;',line ):
		outTop.write( line )
		continue

	if re.search('\[\satomtypes\s\]',line ):
		inAtomTypes = True
		outTop.write( "#include \"./ildnGO.itp\"\n" )
		outTop.write( "#include \"/usr/local/gromacs/share/gromacs/top/amber99sb-ildn.ff/ffbonded.itp\"" )

		continue
	
	if re.search( '\[', line ):
		inAtomTypes = False	

	if not inAtomTypes:
		outTop.write( line )
	
	if re.search( '\[\satoms\s\]', line ):
		break

# In atoms now, so translate the names!
inAtoms = True
for line in top:
	if re.search( '^\s*;',line):
		outTop.write( line )
		continue
	if line.strip() == '':
		inAtoms = False

	if inAtoms:
		atom = line.split()[1]
		res = line.split()[3]
		try: newAtom = RenameDict[ res ][ atom ] # Rename is a dictionary of residue dictionaries...
		except: 
			try: newAtom = CatchAllDict[ atom ]
			except: 
				print "No atom found for %s..." % atom
				exit()
		outTop.write( line.split( atom )[0] + newAtom + atom.join( line.split( atom )[1:] ) )
		continue
	else:
		outTop.write( line )
		if re.search( 'angles', line ):
			break


# Now  in angles, so delete the last two entries so that it's set to the default numbers
for line in top:
	mObj = re.search( '^(\s+\d+\s+\d+\s+\d+\s+\d+\s)', line )
	if not mObj:
		outTop.write( line )
		if re.search( 'dihedrals', line ):
			break
	else:
		outTop.write( mObj.group(1) + ' \n' )
	
# Now in dihedrals, so delete last few entries, but rewrite 1 -> 9 ( proper ) and 2 -> 4 ( improper )
dihTypeDict = { '1': '9', '2': '4' }
for line in top:
	mObj = re.search( '^(\s+\d+\s+\d+\s+\d+\s+\d+\s+)', line )
	if not mObj:
		outTop.write( line )
		if line.strip() == '':
			break
			# this means that you are at the end of the dihedrals section and therefore the end of the topology
	else:
		outTop.write( mObj.group(1) + dihTypeDict[ line.split()[4] ] + ' \n' )

for line in top:
	outTop.write( line )
