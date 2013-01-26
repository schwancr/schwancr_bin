#!/usr/bin/env python

from msmbuilder import Project
import os
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-p',dest='projFN',default='./ProjectInfo.h5')
parser.add_option('-o',dest='outFN',default='./ProjectInfo.h5')

options, args = parser.parse_args()

p = Project.Project.LoadFromHDF( options.projFN )

p['TrajFilePath'] = os.path.abspath( p['TrajFilePath'] )

p['ConfFilename'] = os.path.abspath( p['ConfFilename'] )
os.system('mv ProjectInfo.h5 ProjectInfo.h5.OLD')
p.SaveToHDF(options.outFN)
