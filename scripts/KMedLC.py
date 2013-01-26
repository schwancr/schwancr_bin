#!/home/schwancr/Installed/epd/bin/python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t',dest='Tdata_dir',help='Directory to look for lh5\'s of XYZData over time')
parser.add_option('-u',type=int, default=1,dest='stride', help='Stride to sample by')
parser.add_option('-r',dest='cutoff',type=float, default=-1,help='Cutoff to quit clustering')
parser.add_option('-i',dest='aind_FN',default='../AtomIndices.dat', help='Atom index file to use in calculating RMSD')
parser.add_option('-k',dest='numGens',type=int,default=1000000,help='Maximum number of Generators')
parser.add_option('-m',dest='k_med_iters',type=int,default=10,help='Number of local K-Medoid iterations')
parser.add_option('--cr',dest='coef_rmsd',type=float, default=1, help='Coefficient for rmsd [ 1 ]' )
parser.add_option('--cq',dest='coef_qnorm',type=float, default=1, help='Coefficeint for Q-norm [ 1 ]')
parser.add_option('--ri',dest='rind_FN',help='Residue indices to use in calculating contacts.')
options, args = parser.parse_args()

import numpy as np
from Emsmbuilder import clustering, metrics, Trajectory, Serializer
from pyschwancr import dataIO
import os
import re
Aind = np.loadtxt(options.aind_FN,int)
if options.rind_FN:
   Rind = np.loadtxt(options.rind_FN,int) - 1
else:
   Rind = 'all'

print "Loading XYZData Trajectories ..."
TrajList = dataIO.getTrajList( options.Tdata_dir )
Traj = clustering.concatenate_trajectories( [ Trajectory.Trajectory.LoadFromLHDF( fn )[::options.stride] for fn in TrajList ] )

print "Preparing metrics"
Coefficients = [ options.coef_rmsd, options.coef_qnorm ]
RMSD = metrics.RMSD( atomindices=Aind )
QNorm = metrics.BooleanContact( contacts=Rind, cutoff=0.6 )
LC = metrics.Hybrid( [ RMSD, QNorm ], Coefficients  )
print "Prepared metrics. Now clustering."

KMed = clustering.HybridKMedoids( LC, Traj, options.numGens, distance_cutoff = options.cutoff, local_num_iters= options.k_med_iters )
print "Done clustering. Saving data and exiting."
GenInd = KMed.generator_indices
Ass = KMed.assignments
Dist = KMed.distances

print GenInd

Serializer.SaveData( 'Assignments.h5', Ass )
Serializer.SaveData( 'Assignments.h5.Dists', Dist )
GenTraj = Traj
GenTraj['XYZList'] = GenTraj['XYZList'][GenInd]
GenTraj.SaveToLHDF('gens.lh5')


