#!/usr/bin/env python

from msmbuilder import arglib
from msmbuilder import Project, io, msm_analysis
from scipy.io import mmread
import numpy as np

MaxInt = np.iinfo(np.int32).max

def run( tProb, start_state, steps, project, out_pdb, out_xtc ):

    state_traj = msm_analysis.sample( tProb, start_state, steps )
    print "Sampled tProb."
    state_sizes = np.bincount( assignments[np.where(assignments!=-1)] ) # size of each state.
    which_ind_traj = [ np.random.randint( state_sizes[i] ) for i in state_traj ] # Random integer in each state in the traj

    uniq_states = np.unique( state_traj )

    state_lookup = dict(zip(uniq_states, np.arange(uniq_states.shape[0])))

    print "Translating to trajectory, frame pairs..."
    which_states = [ np.array( np.where( assignments == i ) ).T for i in uniq_states ]
        # The (traj,frame) list for each state
    which_traj = np.array([which_states[state_lookup[uniq_states[state_lookup[state]]]][i] 
                           for state, i in zip(state_traj, which_ind_traj)])
        # Grab a random (traj,frame) for each state visited in the trajectory
    
    print "Loading frames from the Trajectory."

    traj = project.load_frame(which_traj[0][0], which_traj[0][1])

    for i in range(1, len(which_traj)):
        traj += project.load_frame(which_traj[i][0], which_traj[i][1])
        print i
    
    traj[0].save_to_pdb(out_pdb)
    traj.save_to_xtc(out_xtc)

    print "Saved output to %s and %s" % (out_pdb, out_xtc)

    return

if __name__ == '__main__':
    
    parser = arglib.ArgumentParser(description='''
        Sample the transition probability matrix and output
        an xtc file, as well as the pdb of the first state.
        (So you can load it in VMD).''')
    parser.add_argument('start_state',help='Starting state. Pass random for a random state',
        default='random')
    parser.add_argument('tProb')
    parser.add_argument('output',help='Output filename to save xtc and pdb')
    parser.add_argument('project')
    parser.add_argument('assignments')
    parser.add_argument('steps',type=int,help='Number of steps to sample',default=100)
    
    args = parser.parse_args()

    if args.start_state.lower() == 'random':
        start_state = None # msm_analysis will generate the random int
    else:
        start_state = int( args.start_state )

    if args.output[-4:] in ['.xtc','.pdb']:
        out_pdb = args.output[:-4] + '.pdb'
        out_xtc = args.output[:-4] + '.xtc'
    else:
        out_pdb = args.output + '.pdb'
        out_xtc = args.output + '.xtc' 

    project = Project.load_from( args.project )
    assignments = io.loadh(args.assignments)
    try: assignments = assignments['arr_0']
    except: assignments = assignments['Data']

    tProb = mmread( args.tProb )

    run( tProb, start_state, args.steps, project, out_pdb, out_xtc )
