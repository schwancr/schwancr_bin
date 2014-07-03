#!/usr/bin/env python

from msmbuilder import io
from msmbuilder import msm_analysis
from scipy.io import mmread
from argparse import ArgumentParser
import os

parser = ArgumentParser()
parser.add_argument('-t', dest='tProb', help='transition matrix', default='./tProb.mtx')
parser.add_argument('-o', dest='output', help='output filename', default='./eigs.h5')
parser.add_argument('-n', dest='num_vecs', help='number of eigenvectors to find.', default=500, type=int)

args = parser.parse_args()

if os.path.exists(args.output):
    raise Exception("path (%s) exists!" % args.output)

tProb = mmread(args.tProb)

eigs = msm_analysis.get_eigenvectors(tProb, args.num_vecs)

io.saveh(args.output, vals=eigs[0], vecs=eigs[1])
