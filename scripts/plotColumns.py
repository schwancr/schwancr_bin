#!/usr/bin/env python -u 

from optparse import OptionParser
from pyschwancr import dataIO
parser = OptionParser()
parser.add_option('-d',dest='data_FN', help = "Data to plot (col 0 on x-axis, col 1 on y-axis)" )
parser.add_option('-y',dest='y_lbl',help="Y-label")
parser.add_option('-x',dest='x_lbl',help="X-label")
parser.add_option('-o',dest='out_FN',help="Output filename")
options, args = parser.parse_args()


dat = dataIO.readData( options.data_FN )

dataIO.plotData( dat[:,0], dat[:,1], options.x_lbl, options.y_lbl, outName = options.out_FN )

