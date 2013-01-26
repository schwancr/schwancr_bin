#!/usr/bin/env python
from msmbuilder import arglib, io, tICA, Project

parser = arglib.ArgumentParser( get_basic_metric=True )

parser.add_argument( 'project' )
parser.add_argument( 'output' )
parser.add_argument( 'stride', type=int, default=1 )

args, metric = parser.parse_args()

project = Project.load_from( args.project )

arglib.die_if_path_exists( args.output )

stride = int( args.stride )

cov_mat = tICA.CovarianceMatrix( 0 )

for i in xrange( project.n_trajs ):

    print "Working on Trajectory %d"
    ptraj = metric.prepare_trajectory( project.load_traj( i, stride=stride ) )
    if i == 0:
        cov_mat.set_size( ptraj.shape[1] )

    cov_mat.train( ptraj )

print "Saving matrix to %s" % args.output
io.saveh( args.output, covariance_matrix=cov_mat.get_current_estimate() )
