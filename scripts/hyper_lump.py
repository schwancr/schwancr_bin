#!/usr/bin/env python

from argparse import ArgumentParser
from scipy.io import mmread
from schwancrtools import add_lumping
import numpy as np
from msmbuilder import Serializer

parser = ArgumentParser()
parser.add_argument('-t',dest='tProb',help='tProb.mtx from msmbuilder')
parser.add_argument('-k',dest='num_vecs',help='Number of eigenvectors to use')
parser.add_argument('-a',dest='ass_fn',help='Assignments Filename')
args = parser.parse_args()

t = mmread( args.tProb )
k = int( args.num_vecs )
ass = Serializer.LoadData( args.ass_fn )

mapping = add_lumping.hypercube_lump( t, k )

new_ass = np.ones( ass.shape ) * -1

for i,row in enumerate(ass):
    new_ass[i][ np.where( row != -1 ) ] = mapping[ row[ np.where( row != -1 ) ] ]

print np.unique(new_ass)

Serializer.SaveData('HyperMacro.h5', new_ass.astype(int))

