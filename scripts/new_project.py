#!/usr/bin/env python
from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument('-p',dest='proj_FN',help='ProjectInfo.h5 from msmbuilder')
parser.add_argument('-w',dest='write_dir',help='Output directory to create a new project')
parser.add_argument('--min-length',dest='min_length',default=0,type=int,help='Minimum length to keep a trajectory in the new project (in frames)')
parser.add_argument('--trim-first',dest='trim_first',default=0,type=int,help='Trim each trajectory from the beginning by N frames. NOTE: This is applied before min_length is.')

args = parser.parse_args()

if args.min_length < 0:
   args.min_length = 0

if ( args.min_length == 0 ) and ( args.trim_first == 0 ):
   print "Must have a positive minimum length or trim frames for this to be useful"
   exit()

from msmbuilder import Project, Serializer, io
import numpy as np
import tables
import os, re, sys

if os.path.exists( args.write_dir ):
   print "Directory already exists... Delete it before proceeding"
   exit()

os.mkdir( args.write_dir )

#Proj = Project.LoadFromHDF( args.proj_FN )
#Proj = Project.load_from_hdf( args.proj_FN )
Proj = Project.load_from( args.proj_FN )

os.system( 'cp %s %s' % ( Proj.conf_filename, os.path.join( args.write_dir, Proj.conf_filename.split('/')[-1]) ) )

traj_lens = Proj.traj_lengths

not_too_short_inds = np.where( traj_lens >= ( args.min_length + args.trim_first ) )[0]

os.mkdir( os.path.join( args.write_dir, 'Trajectories' ) )
print "Will limit this project to %d trajectories." % len( not_too_short_inds )
for i in xrange( len( not_too_short_inds ) ):
   print "Copying trajectory %d -> %d (length=%d)" % ( not_too_short_inds[i], i, Proj.traj_lengths[ not_too_short_inds[i] ] - args.trim_first )
   trj0 = tables.openFile( Proj.traj_filename( not_too_short_inds[i] ) )
   trj1 = tables.openFile( os.path.abspath( os.path.join( args.write_dir, 'Trajectories', '%s%d%s'% ('trj',i, '.lh5' ) ) ), 'w' )
   #os.system( 'ln -s %s %s' % ( trj0, trj1 ) )
   #os.symlink( trj0, trj1 )
   
   for n0 in trj0.iterNodes('/'):
      if n0.name != 'XYZList':
         trj0.copyNode( where='/', name=n0.name, newparent=trj1.root )
      else:
         temp_ary = n0[ args.trim_first : ]
         io.saveh( trj1, XYZList=temp_ary )

   trj0.close()
   trj1.close()

new_records = {'conf_filename': Proj.conf_filename.split('/')[-1], 
           'traj_lengths': Proj.traj_lengths[ not_too_short_inds ] - args.trim_first,
           'traj_paths': Proj._traj_paths[ : len( not_too_short_inds ) ], # This works because they're named relatively and they are re-numbered
           'traj_converted_from': Proj._traj_converted_from[ not_too_short_inds ],
           'traj_errors': Proj._traj_errors[ not_too_short_inds ] }
new_proj_dir = args.write_dir
# Copy the trajectories
New_Proj = Project( new_records, project_dir = new_proj_dir )
New_Proj.save( os.path.join( args.write_dir, 'ProjectInfo.yaml' ) )
