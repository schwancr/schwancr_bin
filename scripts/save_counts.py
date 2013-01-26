#!/usr/bin/env python
from msmbuilder import io, MSMLib
from scipy.io import mmwrite
import sys
from msmbuilder import arglib

parser = arglib.ArgumentParser()
parser.add_argument('assignments')
parser.add_argument('lagtime', type=int, default=1)
parser.add_argument('sliding_window', default=False, action='store_true')
parser.add_argument('output', default='tCounts.raw.mtx')
args = parser.parse_args()

try: ass = io.loadh(args.assignments)['Data']
except: ass = io.loadh(args.assignments)['arr_0']

C = MSMLib.get_count_matrix_from_assignments( ass, lag_time=args.lagtime, sliding_window=args.sliding_window )

print C.sum()
mmwrite(args.output, C)
