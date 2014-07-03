
import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from msmbuilder import io
from msmbuilder.reduce import tICA
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f1', dest='input1', help='input file #1')
parser.add_argument('-f2', dest='input2', help='input file #2')
parser.add_argument('-n', dest='num_vecs', type=int, help='number of vectors to compare')
parser.add_argument('-o', dest='output', help='output filename')
parser.add_argument('-t1', dest='lag1', type=int, help='lag for file #1')
parser.add_argument('-t2', dest='lag2', type=int, help='lag for file #2')

args = parser.parse_args()

tic1 = io.loadh(args.input1)
tic2 = io.loadh(args.input2)

ind1 = np.argsort(tic1['vals'])[::-1]
ind2 = np.argsort(tic2['vals'])[::-1]

tic1['vals'] = tic1['vals'][ind1]
tic1['vecs'] = tic1['vecs'][:, ind1]

tic2['vals'] = tic2['vals'][ind2]
tic2['vecs'] = tic2['vecs'][:, ind2]

sim_mat = np.zeros((args.num_vecs, args.num_vecs))

for i in xrange(args.num_vecs):
    for j in xrange(args.num_vecs):
        v1i = tic1['vecs'][:, i]
        v2j = tic2['vecs'][:, j]

        temp = v1i.dot(v2j)
        temp /= np.sqrt(v1i.dot(v1i))
        temp /= np.sqrt(v2j.dot(v2j))

        sim_mat[i, j] = temp
 
sim_mat = np.abs(sim_mat)

times1 = - args.lag1 / np.log(tic1['vals'][:args.num_vecs])
times2 = - args.lag2 / np.log(tic2['vals'][:args.num_vecs])
 
figure(figsize=(9, 3))  
ax = subplot(131)
im = ax.matshow(sim_mat, cmap=matplotlib.cm.hot_r)
ylabel('lag=%d' % args.lag1)
xlabel('lag=%d' % args.lag2)
yticks([])
xticks([])

colorbar(im)

subplot(133)

scatter(times2, times1)
ylabel('lag=%d' % args.lag1)
xlabel('lag=%d' % args.lag2)
yscale('log')
xscale('log')

print "saved output to %s" % args.output
savefig(args.output)
