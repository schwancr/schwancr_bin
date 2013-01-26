#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-c', dest='gro_FN', help='Gro file to calculate the contacts for')
parser.add_option('-o',dest='out_FN',help='Output file to write. Really just a basename. This script will create a matrix as well as a list of contacts')
parser.add_option('--co',dest='cutoff',default=0.4, type=float, help='Cutoff to use. Will output 0, 1, 2 for contacting = 0, contacting <= cutoff, contacting > cutoff')
parser.add_option('-m',dest='min_dist',default=4,type=int, help='Minumum distance along a chain to consider something as contacting')

options, args = parser.parse_args()
 
from numpy import *
import os, sys, re
 
# First use g_mdmat to generate the matrix:

mdmatCmd = "echo \"C-alpha\" | g_mdmat -f %s -s %s -mean tmp.xpm -t %f -nlevels 3 &> /dev/null" % ( options.gro_FN, options.gro_FN, options.cutoff * 2 ) # Need to multiply by 2 to get g_mdmat to cooperate

os.system( mdmatCmd )

transDict = { 'A' : 0, 'B' : 1, 'C' : 2 }

xpm = open( 'tmp.xpm' , 'r' )
Rows = []
for line in xpm:
	m = re.search( '^\"([ABC]+)\"', line )
	n = re.search( '^\"(\d+)\s+(\d+)', line )
	if m:
		row = m.group(1)
		Rows.append( row )
	elif n:
		dim = int( n.group(1) )
		if dim != int( n.group(2) ):
			print "Unequal dimensions!!!"
			exit()
	else:
		continue

Mat = zeros( ( dim, dim ) )

for i in range( dim ):
	for j in range( i, dim ):
		Mat[ i, j ] = transDict[ Rows[i][j] ]
		Mat[ j, i ] = Mat[ i, j ]

savetxt( options.out_FN +'_ContactMap.dat', Mat )

contacting = array(  where( Mat < 2 ) ).T
contacting += 1
contactFile = open( options.out_FN + '.contacts', 'w' )

for a,b in contacting:
	if ( b - a ) > options.min_dist:
		contactFile.write( '1 %d 1 %d \n' % ( a, b ) )
		
contactFile.close()
xpm.close()

