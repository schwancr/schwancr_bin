#!/usr/bin/env python
from msmbuilder import arglib
from msmbuilder import msm_analysis, MSMLib
from scipy.io import mmread
import scipy.sparse
import multiprocessing as mp
import numpy as np

def get_eigenvalues( count_matrix ):


    bad_states = np.array(np.where( count_matrix.sum(axis=1) == 0 )[0]).flatten()

    i_ary = count_matrix.nonzero()[0]
    j_ary = count_matrix.nonzero()[1]

    i_ary = np.concatenate( (i_ary, bad_states) )
    j_ary = np.concatenate( (j_ary, bad_states) )
    new_data = np.concatenate( (count_matrix.data, np.ones(len(bad_states))) )

    print i_ary.shape, count_matrix.data.shape, new_data.shape, len(bad_states)

    count_matrix = scipy.sparse.csr_matrix( (new_data, (i_ary, j_ary)) )

    #count_matrix = count_matrix.tolil()
    #count_matrix[(bad_states, bad_states)] = 1
    #count_matrix = count_matrix.tocsr()

    print count_matrix.data.shape, count_matrix.nonzero()[0].shape
    #NZ = np.array(count_matrix.nonzero()).T

    #keep_ind = []
    #for i in xrange(len(NZ)):
    #    if NZ[i][0] in bad_states or NZ[i][1] in bad_states:
    #        pass
    #    else:
    #        keep_ind.append(i)
    #keep_ind = np.array(keep_ind)

    #N = NZ.max()+1

    #count_matrix = scipy.sparse.csr_matrix( (np.array(count_matrix.data)[keep_ind], NZ[keep_ind].T), shape=(N,N), copy=True )

    try:
        t_matrix = MSMLib.build_msm(count_matrix, symmetrize=args.symmetrize)[1]
    except:
        return None
    
    vals = msm_analysis.get_eigenvectors(t_matrix, args.num_vals, epsilon=1)[0]
    vals.sort()

    return vals[::-1]

parser = arglib.ArgumentParser()
parser.add_argument('counts', default='tCounts.raw.mtx', 
                    help='Transition counts before any trickery.')
parser.add_argument('num_vals', default=10, type=int, 
                    help='Number of eigenvalues to analyze')
parser.add_argument('output', default='bootstrap_err.dat')
parser.add_argument('symmetrize', help="""Method by which to estimate a
    symmetric counts matrix. Symmetrization ensures reversibility, but may skew
    dynamics. We recommend maximum likelihood estimation (MLE) when tractable,
    else try Transpose. It is strongly recommended you read the documentation
    surrounding this choice.""", default='MLE',
    choices=['MLE', 'Transpose', 'None'])
#parser.add_argument('lag_time', type=int)
parser.add_argument('num_samples', default=1000, help="""Number of times to
    resample the counts matrix.""", type=int)
parser.add_argument('procs', default=1, type=int)
#parser.add_argument('num_steps',type=int)

args = parser.parse_args()

raw_counts = mmread(args.counts)

sampled_c_matrices = msm_analysis.resample_c_matrix(raw_counts, 
                                                    num_samples=args.num_samples)

    
#sampled_trajs = [np.array(msm_analysis.sample(raw_counts, 0, args.num_steps)) for i in range(args.num_samples)]
#sampled_c_matrices = [MSMLib.get_counts_from_traj(s, lag_time=args.lag_time, sliding_window=False) for s in sampled_trajs]

pool = mp.Pool(args.procs)

result = pool.map_async(get_eigenvalues, sampled_c_matrices)
result.wait()
eigenvalues = result.get()

eigenvalues = [ v for v in eigenvalues if not v is None ]

eigenvalues = np.array(eigenvalues)

print "We had a total of %d good samples."%len(eigenvalues)
stdevs = eigenvalues.std(axis=0)
means = eigenvalues.mean(axis=0)

np.savetxt(args.output, np.array((means, stdevs)).T)

