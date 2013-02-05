
from msmbuilder import arglib
from msmbuilder import io
from msmbuilder import utils
from msmbuilder import Project
from msmbuilder import Trajectory
import numpy as np
import os, sys, re

parser = arglib.ArgumentParser()
parser.add_argument('traj_dir', help='Directory to find trajectory files.')
parser.add_argument('conf_fn', help='Conformation filename that has the same atom names and residue IDs, etc. as the trajectories.')
parser.add_argument('output', default='./ProjectInfo.yaml', help='Output filename [ ./ProjectInfo.yaml ]')

args = parser.parse_args()

traj_list = [ os.path.join(args.traj_dir, fn) for fn in os.listdir(args.traj_dir)]

traj_list.sort(key=utils.keynat) # = list.sort(traj_list, key=utils.keynat)
print traj_list

traj_lens = []

for i in xrange(len(traj_list)):
    
    print i
    shape = Trajectory.load_from_hdf(traj_list[i], JustInspect=True)
    traj_lens.append(shape[0])

records = { 'conf_filename' : args.conf_fn,
            'traj_lengths' : traj_lens,
            'traj_paths' : traj_list,
            'traj_converted_from' : [[] for fn in traj_list],
            'traj_errors' : [None for fn in traj_list]
          }

project = Project(records=records)

project.save(args.output)
