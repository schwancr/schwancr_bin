#!/usr/bin/env python

from msmbuilder import Trajectory
import numpy as np
import os, re, sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-t',dest='traj',help='Trajectory to subsample')
parser.add_argument('-o',dest='out',help='Output filename')
parser.add_argument('-N',dest='num_frames',type=int,help='Number of frames to keep at the beginning')

options = parser.parse_args()

trj = Trajectory.LoadTrajectoryFile( options.traj )

if len( trj ) <= options.num_frames:
   print "Trajectory's length (%d) is too short for your input (%d)." % ( len(trj), options.num_frames )
   exit()

trj = trj[:options.num_frames]

trj.SaveToLHDF( options.out )
