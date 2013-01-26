#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-d', dest='dirs',action='append',help='Directory to find a ProjectInfo.h5 and Trajectories directory. Specify more than one to combine.')
parser.add_option('-w', dest='writeDir',default='Combined',help='Directory to write ProjectInfo.h5 and Trajectories')

options, args = parser.parse_args()

import os, shutil
from msmbuilder import Serializer, Project
import numpy as np

# First read in the first project, we will merely add to this template

Proj = Serializer.Serializer.LoadFromHDF( os.path.join( options.dirs[0], 'ProjectInfo.h5' ) )

NumTrajs = [ Proj['NumTrajs'] ]
print "Combining project info files"
for dir in options.dirs[1:]:
   ProjNew = Serializer.Serializer.LoadFromHDF( os.path.join( dir, 'ProjectInfo.h5' ) )
   Proj['TrajLengths'] = np.concatenate( ( Proj['TrajLengths'], ProjNew['TrajLengths'] ) )
   Proj['RunList'] = np.concatenate( ( Proj['RunList'], ProjNew['RunList'] ) )
   Proj['NumGensList'] = np.concatenate( ( Proj['NumGensList'], ProjNew['NumGensList'] ) )
#   Proj['TrajFileBaseName'] Same for all projects...
#   Proj['TrajFilePath'] Same for all projects ... since it's a relative directory
#   Proj['ConfFilename'] Same for all projects ... since it's a relative directroy... This could cause problems
#   Proj['TrajFileType'] Same for all projects... presumably... lh5
   Proj['CloneList'] = np.concatenate( ( Proj['CloneList'], ProjNew['CloneList'] ) )
   Proj['NumTrajs'] = Proj['NumTrajs'] + ProjNew['NumTrajs']
   NumTrajs.append( ProjNew['NumTrajs'] )


if not os.path.isdir( options.writeDir ):
   os.mkdir( options.writeDir )
  
NewTrajDir = os.path.join( options.writeDir, 'Trajectories' ) 
if os.path.isdir( NewTrajDir ):
   print "Trajectories directory exists! Remove this and run the script again..."
   exit()

# Directory exists, so now copy the trajectories.
print "Copying data to the new Trajectories directory:"

readmeFN = open( os.path.join( options.writeDir, 'README' ),'w')
readmeFN.write( 'This directory was created by combining projects in (in the indicated order):\n')
for dir in options.dirs:
   readmeFN.write( os.path.abspath(dir) + '\n' )
readmeFN.close()

#shutil.copytree( os.path.join( options.dirs[0], 'Trajectories' ), NewTrajDir )
# shutil is SLOW.... ^^^^ so use os.system

os.system( 'cp -r %s %s' % ( os.path.join( options.dirs[0], 'Trajectories' ), NewTrajDir ) )

sumInd = 0
for dirInd,dir in enumerate( options.dirs ):
   if dirInd == 0:
      sumInd += NumTrajs[ dirInd ]
      continue

   for trjIndTemp in range( NumTrajs[ dirInd ] ):
      #shutil.copy( os.path.join(dir,'Trajectories','trj%d.lh5' % trjIndTemp ), os.path.join( NewTrajDir, 'trj%d.lh5' % ( trjIndTemp + sumInd ) ) )
      # shutil is SLOW so use os.system
      os.system( 'cp %s %s' % ( os.path.join(dir,'Trajectories','trj%d.lh5' % trjIndTemp ), os.path.join( NewTrajDir, 'trj%d.lh5' % ( trjIndTemp + sumInd ) ) ) )
      #print os.path.join(dir,'Trajectories','trj%d.lh5' % trjIndTemp ), os.path.join( NewTrajDir, 'trj%d.lh5' % ( trjIndTemp + sumInd ) )
      print "Project %d ..." % dirInd, "Trajectory %d\r" % trjIndTemp,
   sumInd += NumTrajs[ dirInd ]

print "\n"
# Proj now represents the project for the new combined directory.
Proj.SaveToHDF( os.path.join( options.writeDir, 'ProjectInfo.h5' ) )

# Move the trajectories.

os.chdir( options.writeDir )
try: p = Project.Project.LoadFromHDF('ProjectInfo.h5')
except Exception as ERR:
   print "Something went wrong in writing the ProjectInfo.h5 file. Here's the error when it was loaded:"
   print ERR
   print "You should fix this issue yourself."



