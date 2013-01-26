#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-a','--assignment-file',dest='assFN',help='Assignments.h5 containing the assignments for your MSM')
parser.add_option('-p','--project-file',dest='projFN',help='Project filename, to use in reshaping the raw data')
parser.add_option('-d','--data-file',dest='dataFN',help='Numpy array containing the native contacts for all snapshots')
parser.add_option('-o','--output-file',dest='outFN',default='AssignmentNCs.dat',help='Output to write to (flat text)')
parser.add_option('-P','--procs',dest='procs',type='int',default=1,help='Number of threads/processors to use.')
options,args = parser.parse_args()

from msmbuilder import io, Project
from numpy import *
import os
import sys
import multiprocessing as mp
from pyschwancr import msmTools

global count
count = 0
def calcAvg( cluster ):
    global count
    global options
    count += 1
    print "\r%6.2f %% Clusters Finished..." % (float(count) * int(options.procs) * 100 / (As.max()+1),),
    ind = where( As == cluster )
    avg = data2d[ind].mean()
    SD = data2d[ind].std()
	
    return (cluster, avg, SD)

def main():
    global data2d
    global As
    # First I need to turn the assignments matrix into a 1D list of assignments
    sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)
    print "Reading in Assignments... from %s " % options.assFN
    As = io.loadh(options.assFN)['arr_0'].astype(int)
    proj = Project.load_from( options.projFN )
    print "Reading in data... from %s " % options.dataFN
    try: 
        f = io.loadh( options.dataFN )
        try:
            data2d = f['arr_0']
        except:
            data2d = f['Data']
    except:
        data = load(options.dataFN)
        data2d = msmTools.reshapeRawData( data, proj )

    print "Calculating averages for:"

    pool = mp.Pool(options.procs)
    clusters = range( As.max() + 1)
    result = pool.map_async(calcAvg,clusters[:])
    result.wait()
    sol = result.get()
    sol = array(sol)
    savetxt(options.outFN, sol)

    return

if __name__=='__main__':
    main()
