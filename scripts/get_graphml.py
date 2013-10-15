
from scipy.io import mmread
import networkx as nx
import sys

f = sys.argv[1]

fout = f[:-4] + '.graphml'

d = mmread(f)

g = nx.from_scipy_sparse_matrix(d)

nx.write_graphml(g, fout)
print 'saved output to %s' % fout
