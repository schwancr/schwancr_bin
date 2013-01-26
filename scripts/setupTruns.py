#!/usr/bin/env python 


from optparse import OptionParser
import sys
import re
import os

parser = OptionParser()
parser.add_option('-p',dest='top_FN',help="Gromacs topology file" )
parser.add_option('-c',dest='gro_FN',help="Gromacs filename" )
parser.add_option('-w',dest='out_Dir',default='./Sims/',help="Directory to place tpr's [ ./Sims/ ]")
parser.add_option('-f',dest='mdp_FN',help="Gromacs mdp file")
parser.add_option('--t0',dest='temp_0',type=float, default=100., help="Low temp [ 100 ]" )
parser.add_option('--tf',dest='temp_f',type=float, default=200., help="High temp [ 200 ]" )
parser.add_option('-n',dest='num_runs',type=int,default=60, help="Number of tpr's to generate [ 60 ]" )
parser.add_option('--zero',dest='zero',action='store_true',default=False,help="Pass this flag to make non-default bond/angle/torsions set to zero rather than throwing an error during grompp")
options, args = parser.parse_args()

print sys.argv[1]

os.putenv('GMX_MAXBACKUP','-1')
init_mdp = open( options.mdp_FN ,'r')

mdpLines = init_mdp.readlines()

workDir = os.getcwd()
top = options.top_FN
gro = options.gro_FN

N = options.num_runs
temp0 = options.temp_0
tempF = options.temp_f

if not os.path.exists( options.out_Dir ):
   os.system( 'mkdir %s' % options.out_Dir )
if not os.path.exists( os.path.join( options.out_Dir, 'Mdps' ) ):
   os.system( 'mkdir %s' % os.path.join( options.out_Dir, 'Mdps' ) )


for i in range(N):

   mdpTemp = open(os.path.join( options.out_Dir, 'Mdps',"mdpTemp"+str(i)+".mdp") ,'w')

   temp = temp0 + (tempF - temp0)*i/float(N-1)

   for line in mdpLines:
      if re.search('ref_t',line):
         outLine = 'ref_t     = '+str(temp)+'\n'
      elif re.search('gen_temp',line):
         outLine = 'gen_temp  = '+str(temp)+'\n'
      else:
         outLine = line

      mdpTemp.write(outLine)

   mdpTemp.close()

   cmd = "grompp -f %s -p %s -c %s -o %s -maxwarn 1" % ( mdpTemp.name, top, gro, os.path.join( options.out_Dir, str(i) ) )
   if options.zero:
      cmd += " -zero"
   os.system( cmd )

