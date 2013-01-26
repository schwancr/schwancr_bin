#!/home/schwancr/Installed/epd/bin/python
# This file is part of MSMBuilder.
#
# Copyright 2011 Stanford University
#
# MSMBuilder is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

import sys
import os
import glob
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t',dest='Rdata_dir',help='Data directory to find XYZData for trajectories')
parser.add_option('-q',dest='Qdata_dir',help='Data directory to find QData for trajectories')
parser.add_option('--cr',dest='coef_rmsd',default=1,type=float, help='Coefficient for RMSD [ 1 ]')
parser.add_option('--cq',dest='coef_qnorm',default=1,type=float, help='Coefficient for Q-norm [ 1 ]')
parser.add_option('--gq',dest='Qgens_FN',default='./Gens_QData.npy',help='Generator\'s native contact array [ Gens_QData.npy ]')
parser.add_option('--gr',dest='Rgens_FN',default='./Gens_XYZData.lh5',help='Generator\'s XYZ trajectory [ Gens_XYZData.lh5 ]')
parser.add_option('-P',dest='procs',default=1,type=int,help='Number of procs used if doing multiple assigning')
parser.add_option('-i',dest='aind',help='Atom indices to use in calculating RMSDs')
parser.add_option('--which-queue',dest='whichQueue',default='long',help='Which queue to submit pbs scripts to [ long ]')
options, args = parser.parse_args()


def run(XYZ_TrajDir, Q_TrajDir, XYZ_GenFN, Q_GenFN, XYZ_Coef, Q_Coef, IndFilename, procs, WhichQueue):
  # Check Output
  if len(glob.glob('./AssignOnPBS*.sh')) != 0:
    print "Error: There are already 'AssignOnPBS*.sh' files present! Exiting."
    sys.exit(1)
  
  # Create a bunch of AssignOnCertaintyXX.sh scripts and submit them
  for i in range(procs):
    argstring = "-t %s -q %s " % ( XYZ_TrajDir, Q_TrajDir )
    argstring += "--gr %s --gq %s " % ( XYZ_GenFN, Q_GenFN )
    argstring += "--cr %f --cq %f " % ( XYZ_Coef, Q_Coef )
    argstring += "-P %d -N %d " % ( procs, i )
    if IndFilename: # If it's not given then MA_LC.py can handle it.
        argstring += "-i %s " % ( IndFilename )
    cmd='MA_LC.py %s >& AssignPart%d.log' % (argstring, i)
    curDir = os.path.abspath(".")
    name="AssignPart%d" % i

    PBS_File="""#!/bin/bash

#PBS -N %s
#PBS -e /home/schwancr/PBS_LOGFILES/AssignOnPBS%d.sh.err
#PBS -o /home/schwancr/PBS_LOGFILES/AssignOnPBS%d.sh.out
#PBS -l nodes=%d:ppn=%d
#PBS -l walltime=%s
#PBS -V

PBS_O_WORKDIR='%s'
export PBS_O_WORKDIR
### ---------------------------------------
### BEGINNING OF EXECUTION
### ---------------------------------------

        echo The master node of this job is `hostname`
        echo The working directory is `echo $PBS_O_WORKDIR`
        echo This job runs on the following nodes:
        echo `cat $PBS_NODEFILE`


### end of information preamble
cd $PBS_O_WORKDIR

# execute commands
%s
""" % (name, i, i, 1, 24, "20:00:00:00", curDir, cmd)

    fn="AssignOnPBS%s.sh" % i
    f=open(fn, 'w')
    f.write(PBS_File)
    f.close()
    print "Wrote, Submitted, Removed:", fn
    os.system("qsub -q %s %s" % (WhichQueue,fn))
#    os.system("rm %s" % fn)

  return


if __name__ == "__main__":
  print """\nAssigns any data not used in original clustering to generators. Does
this in parallel on a PBS cluster, dividing work amongst many nodes (specify
with -P option). NOTE: The maximum available parallelization is one machine per
trajectory. After this script's children on the nodes completes, run
MergeAssign.py, and the final output is:
-- Assignments.h5: a matrix of assignments where each row is a vector
corresponding to a data trajectory. The values of this vector are the cluster
assignments.
-- Assignments.h5.RMSD: Gives the RMSD from the assigned frame to its Generator.
-- Assignments.h5.WhichTrajs: Shows which trajectories were used (if
subsampling, disabled in this wrapper script).\n"""

  print sys.argv
  
  run( options.Rdata_dir, options.Qdata_dir, options.Rgens_FN, options.Qgens_FN, options.coef_rmsd, options.coef_qnorm, options.aind, options.procs, options.whichQueue )
