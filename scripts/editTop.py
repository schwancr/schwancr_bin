#!/usr/bin/env python

import sys, re, os

top = open( sys.argv[1], 'r' )
mult = float( sys.argv[2] )
outTop = open( sys.argv[1][:-4] + 'NEW.top', 'w')

for line in top:
	outTop.write( line )
	if re.search( '\[ pairs \]', line ):
		break
print line
for line in top:
	mObj = re.search( '(\s*\d+\s*\d+\s*\d+\s*)([\d.E+-]+)\s*([\d.E+-]+)', line )
	if mObj:
		nums = [ mult * float( i ) for i in mObj.groups()[1:] ]
		outLine = mObj.groups()[0] + '%.4E  %.4E \n' % tuple( nums )
		outTop.write( outLine )
	else:
		outTop.write( line )


	if re.search( '\[ bonds \]', line ):
		break

for line in top:
	outTop.write( line )

