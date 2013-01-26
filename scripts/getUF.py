#!/usr/bin/env python

from optparse import OptionParser
parser=OptionParser()
parser.add_option('-f','--order-parameter',dest='dataFN',help='Order parameter to determine unolded and folded ensembles (eg. RMSD to native state)')
parser.add_option('-p','--populations-file',dest='popFN',help='Populations.dat generated from BuildMSM.py')
parser.add_option('--fc','--f-column',dest='fcol',default=1,type=int,help='Columnn in order-param array to use. Default = 1 (indexing from 0)')
parser.add_option('-l','--lower-percent',dest='low',type=float,default=5.,help='Percent to define the lower state (Unfolded)')
parser.add_option('-u','--upper-percent',dest='upper',type=float,default=5.,help='Percent to define the upper state (Folded)')
parser.add_option('-r','--reverse-upper-is-native',dest='upper_is_native',default=True,action='store_false',help="This script assumes a high value in the order parameter is \"Folded\" and a low value is \"Unfolded.\" If this is the opposite then use this flag to make it so.")

options, args = parser.parse_args()

from numpy import *

if options.upper_is_native:
	print "Assuming large order parameter is native"
else:
	print "Large order parameter is non-native"

if options.dataFN.split('.')[-1] == 'npy':
	data = load(options.dataFN)
else:
	print "Loading %s as flat text..." % options.dataFN
	try: 
		data = loadtxt(options.dataFN)
	except: 
		print "Could not load %s... Save data as flat text (numpy savetxt) or as a numpy array and try again." % options.dataFN; exit(1)

pops = loadtxt(options.popFN)


if len(data.shape) > 1:
	print "Data has more than one column. Using column %d. If this is not correct then use --fc to use the correct column" % options.fcol
	data = data[:,options.fcol]

print "All data is loaded."

low = options.low / 100.
up = options.upper / 100.

# Sort the data: First zip it with it's indices
data = zip( data, range( len(data) ) )
# Sort:
data.sort()
lowInd = []
upInd = []
data = array( data )
counter = 0
upNotDone = True
lowNotDone = True
lowPop = 0
upPop = 0
while upNotDone or lowNotDone:

	if counter >= len(data):
		print "Something went wrong..."
		exit(1)

	if lowNotDone:
		lowInd.append(data[counter][1])
		lowPop = pops[ lowInd ].sum()
		if lowPop >= low:
			lowNotDone = False
	
	if upNotDone:
		upInd.append(data[-(counter+1)][1])
		upPop = pops[ upInd ].sum()
		if upPop >= up:
			upNotDone = False
	counter+=1

if options.upper_is_native: # Low = U High = F
	lowFN = "U_%.1fpct_%s.dat" % (options.low, '.'.join(options.dataFN.split('/')[-1].split('.')[:-1] ) )
	upFN = "F_%.1fpct_%s.dat" % (options.low, '.'.join(options.dataFN.split('/')[-1].split('.')[:-1] ) )
	f_cut = data[:,0][ upInd ].min()
	u_cut = data[:,0][ lowInd ].max()
else: # Low = F High = U
	lowFN = "F_%.1fpct_%s.dat" % (options.low, '.'.join(options.dataFN.split('/')[-1].split('.')[:-1] ) )
	upFN = "U_%.1fpct_%s.dat" % (options.low, '.'.join(options.dataFN.split('/')[-1].split('.')[:-1] ) )
	f_cut = data[:,0][ lowInd ].max()
	u_cut = data[:,0][ upInd ].min()


savetxt( 'Cutoffs.dat', array([ f_cut, u_cut ]) )
print "Wrote Cutoffs.dat with the cutoffs used in order folded then unfolded"
savetxt(lowFN,lowInd,"%d")
savetxt(upFN,upInd,"%d")
