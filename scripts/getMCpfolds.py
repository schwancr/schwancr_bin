#!/usr/bin/env python -u

from optparse import OptionParser
parser=OptionParser()
parser.add_option('-x', dest='x_FN', help="X-data filename")
parser.add_option('-y', dest='y_FN', help="Y-data filename")
parser.add_option('-p', dest='proj_FN',help='ProjectInfo.h5 from msmbuilder')
parser.add_option('-r', dest='raw_FN',help='Raw filename with trj##_frm##\nFolded = X Unfolded = Y Neither = Z Total = X+Y+Z')
parser.add_option('-f', dest='f_cut',type=float,help='Folded cutoff. NOTE: greater than f_cut is folded.')
parser.add_option('-u', dest='u_cut',type=float,help='Unfolded cutoff. NOTE: less than u_cut is unfolded.')
parser.add_option('-o', dest='outFN',default='MC_vs_Raw_PFolds.dat', help='Output file to write data to [ ./MC_vs_Raw_PFolds.dat ]')
options, args = parser.parse_args()

from pyschwancr import dataIO, MonteCarlo
from numpy import *
from matplotlib.backends import backend_pdf
from matplotlib.pyplot import *
from msmbuilder import Project
import re

xDat = dataIO.readData(options.x_FN)
yDat = dataIO.readData(options.y_FN)
Proj = Project.Project.LoadFromHDF( options.proj_FN )
rawFN = open( options.raw_FN )

rawText = rawFN.read()

rawLines = rawText.split('trj')

rawXY = []
rawPfolds = []

trajLengths = Proj['TrajLengths']


if len( xDat.shape ) > 1:
	xDat = xDat[:,0]
	print "Using first column of the x-data"
if len( yDat.shape ) > 1:
	yDat = yDat[:,0]
	print "Using first column of the y-data"

for line in rawLines[1:]:
	mObj = re.search( '^(\d+)_frm(\d+)\s*Folded = (\d+) Unfolded = (\d+)', line )
	if mObj:
		trj, frm, F, U = [ float( i ) for i in mObj.groups() ]
		ind = trajLengths[:trj].sum() + int( frm )
	#	print line, trj, frm, F, U
		if F+U < 75:
			print "Less than 75 total simulations folded or unfolded..."
			continue
		rawXY.append( array([ xDat[ ind ], yDat[ ind ] ]) )
		rawPfolds.append( F / ( F + U ) )
	else:
		print line

rawXY = array( rawXY )
print "Setting up the mover ...",
Mover = MonteCarlo.MC_Mover( xDat, yDat )

minXary_F = abs( options.f_cut - Mover.xTicks )
minXary_U = abs( options.u_cut - Mover.xTicks )

fold_ind = where( minXary_F == minXary_F.min() )[0]
unfold_ind = where( minXary_U == minXary_U.min() )[0]
def isF( r ):
	if r[0] > fold_ind:
		return True
	else:
		return False

def isU( r ):
	if r[0] < unfold_ind:
		return True
	else:
		return False

def Either( r ):
	if isF( r ):
		return True
	elif isU( r ):
		return True
	else:
		return False

#r0 = array( [ 0.74025974,  0.8125    ] )

#fold = 0
#tot = 0
#for i in range( 100 ):
	
#	check = MonteCarlo.RandomWalk( Mover, r0, N = 1000, untilFcn = Either )
#	if check:
#		if isF( Mover.positions[-1] ):
#			fold += 1
#		tot += 1
	
#print fold, tot
#exit()

MC_Pfolds = []
for pair in rawXY:
	xyInd = Mover.getXYindices( pair )

	fold = 0
	total = 0	
	for i in range( 100 ):
		check = MonteCarlo.RandomWalk( Mover, xyInd, N=10000, untilFcn=Either )

		if check:
			if isF( Mover.positions[-1] ):
				fold += 1
			total += 1
	print pair, fold, total
	MC_Pfolds.append( float( fold ) / float( total ) )

outAry = array( zip( rawXY[:,0], rawXY[:,1], rawPfolds, MC_Pfolds ) )
print "Saving data as rawX, rawY, rawPfolds, MonteCarloPfolds to %s" % options.outFN
savetxt( options.outFN, array( outAry ) )

#imshow( log10( Mover.Pot.T ), extent=[0,1,0,1], origin='lower')

#dat = Mover.getXYvalues()
#plot( dat[:,0], dat[:,1], color='white' )
#plot( dat[0,0], dat[0,1], 'go')
#plot( dat[-1,0], dat[-1,1], 'rs')

#savefig('randomWalk.pdf')

