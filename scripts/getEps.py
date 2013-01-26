#!/usr/bin/env python

import sys
from numpy import *
import os
import re

def eps( V, W ):
	m = 12.
	n = 10.
	return ( ( V / 4.)**(m/n) * 4./W )**( 1./(m/n-1.) )

f = open( sys.argv[1], 'r')

lines = f.readlines()

for index, line in enumerate(lines):
	if re.search( '\[\spairs\s\]', line):
		break
print index
# At the pair section
Vs = []
Ws = []

for line in lines[index+1:]:
	if len( line.split() ) != 5:
		print "Done: %s" % line
		break
	temp = line.split()
	Vs.append( float( temp[-2] ) )
	Ws.append( float( temp[-1] ) )
	#print temp
epsS = []

for v,w in zip( Vs, Ws ):
	epsS.append( eps( v,w) )
	
savetxt( 'Epsilons_%s.dat' % sys.argv[1][:-4], epsS )
	
print mean( epsS ), "+/-", std( epsS )


