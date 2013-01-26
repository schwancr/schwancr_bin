#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-a','--assignment-file',dest='assFN',help='Assignments.h5 containing the assignments for your MSM')
parser.add_option('-d','--data-file',dest='dataFN',help='Numpy array containing the native contacts for all snapshots')
parser.add_option('-o','--output-file',dest='outFN',default='AssignmentNCs.dat',help='Output to write to (flat text)')
parser.add_option('-P','--procs',dest='procs',type='int',default=1,help='Number of threads/processors to use.')
options,args = parser.parse_args()

from msmbuilder import Serializer
from numpy import *
import os
import sys
import multiprocessing as mp

global count
count = 0
def calcAvg( cluster ):
	global As1D
	global data
	global count
	global options
	count += 1
	print "\r%6.2f %% Clusters Finished..." % (float(count) * int(options.procs) * 100 / As1D.max(),),
	ind = where( As1D == cluster )
	avg = data[ind].mean()
	SD = data[ind].std()
	
	return (cluster, avg, SD)

def main():
	global As1D
	global data
	# First I need to turn the assignments matrix into a 1D list of assignments
	sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)
	As = Serializer.LoadData(options.assFN).astype(int)
	data = load(options.dataFN)

	print "Reading in Assignments... from %s " % options.assFN
	As1D = As[ where( As >= 0 )]
#	As1D = mp.Array('i',As1D)
	#print As1D.shape
	print "Reading in data... from %s " % options.dataFN
	data = load(options.dataFN)
	#print data.shape
	if len(As1D) != len(data):
		print "Assignments (%d) and data (%d) are not the same length!\n \t Exiting..." % (len(As1D), len(data) )
		exit()

	if len(data.shape) > 1:
		print "Data contains more than one column! Using the first column..."
		data = data[:,0]
#	data = mp.Array('d',data)
	print "Calculating averages for:"

	
	pool = mp.Pool(options.procs)
	clusters = range( As1D.max() + 1)
	result = pool.map_async(calcAvg,clusters[:])
	result.wait()
	sol = result.get()
	sol = array(sol)
	savetxt(options.outFN, sol)

	return

if __name__=='__main__':
	main()
