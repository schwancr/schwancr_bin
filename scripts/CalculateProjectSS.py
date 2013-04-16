from msmbuilder import arglib
import sys, re, os
from msmbuilder import Project, Serializer, Trajectory
import numpy as np
import multiprocessing as mp
import subprocess
import warnings
from scipy import array, chararray, sign

#TEMPLOC='./tmp'
#if os.path.exists( TEMPLOC ):
#    warnings.warn('Using %s for temporary files.'%TEMPLOC)
#else:
#    os.mkdir( TEMPLOC )

DSSPBINARY = None # global variable which holds the name of the dssp binary

DEVNULL=open(os.devnull,'w')

CONF=None # Empty container to store the atom names and stuff the PDB writer need

def WritePDB( atom_nums, atoms, res_names, res_nums, xyz, chain ):

    F=''
    for i in range(len(atom_nums)):
        line=chararray(80)
        line[:]=' '
        line[0:4]=array(list("ATOM"))
        line=array(line,'str')
        line[6:11]=array(list(str(atom_nums[i]%100000).rjust(5)))
        #Molprobity is picky about atom name centering
        if len(str(atoms[i]))==3:
            line[12:16]=array(list(str(atoms[i]).rjust(4)))
        elif len(str(atoms[i]))==2:
            line[12:16]=array(list(" "+str(atoms[i])+" "))
        elif len(str(atoms[i]))==1:
            line[12:16]=array(list(" "+str(atoms[i])+"  "))
        else:
            line[12:16]=array(list(str(atoms[i]).center(4)))
        if len(str(res_names[i]))==3:
            line[17:20]=array(list(str(res_names[i])))
        else:
            line[17:21]=array(list(str(res_names[i]).ljust(4)))

        line[21]=str(chain[i]).rjust(1)
        line[22:26]=array(list(str(res_nums[i]).rjust(4)))

        x=10*xyz[i][0]
        y=10*xyz[i][1]
        z=10*xyz[i][2]
        sx=sign(x)
        sy=sign(y)
        sz=sign(z)

        line[30:38]=array(list(("%8.3f"%(x))))
        line[38:46]=array(list(("%8.3f"%(y))))
        line[46:54]=array(list(("%8.3f"%(z))))


        if atom_nums[i]!=-1:
            F += ('%s\n'%line.tostring())
    F += "ENDMDL\n"

    return F

def analyze_conf( xyzlist ):
    #temp_fn = os.path.join( TEMPLOC, str( np.random.randint( 1E9 ) ) + '.pdb' )
    #conf.SaveToPDB( temp_fn )

    conf_str = WritePDB( CONF['AtomID'], CONF['AtomNames'], CONF['ResidueNames'], CONF['ResidueID'], xyzlist, CONF['ChainID'] )

    dssp_output = subprocess.check_output( 'echo "%s" | dssp --' % conf_str, shell=True, stderr=DEVNULL)

    dssp_str = ''.join([ o[16] for o in dssp_output.split('#')[-1].split('\n')[1:] if len(o)>1]).replace(' ','_')

    #os.remove( temp_fn )

    return dssp_str

def run( project, output, num_procs=1, chunk_size=50000 ):

    pool = mp.Pool( num_procs )

    dssp_assignments = []

    for i in xrange( project['NumTrajs'] ):
        traj_dssp_assignments = []
        N = project['TrajLengths'][i]
        j = 0
        for trj_chunk in Trajectory.EnumChunksFromLHDF( project.GetTrajFilename( i ), ChunkSize=chunk_size ):
            result = pool.map_async( analyze_conf, trj_chunk['XYZList'] )
            result.wait()

            traj_dssp_assignments.extend( result.get() )

            j+=len(trj_chunk)
            print "Trajectory %d: %d / %d" % (i, j, N)
        dssp_assignments.append( traj_dssp_assignments )

   # conf['XYZList'] = np.array([ project.ReadFrame(0,0) ])

    dssp_assignments = np.array( dssp_assignments )
    np.save( output, dssp_assignments )
    #Serializer.SaveData( output, dssp_assignments )
    DEVNULL.close()

if __name__ == '__main__':
    parser = arglib.ArgumentParser(
        description='''Calculate the secondary structure assignments
            for all residues in all trajectories according to DSSP.
            NOTE: you must have the dssp binary installed in order 
            for this to work. Pass the name (in the path) or the 
            absolute path with the dssp_binary option''' )

    parser.add_argument('dssp_binary', help='Binary from http://swift.cmbi.ru.nl/gv/dssp/',
        default='dssp')
    parser.add_argument('project')
    parser.add_argument('output',help='Output filename (.h5)',default='DSSP.h5') 
    parser.add_argument('num_procs',help='Number of processes to run.', default=1, type=int)

    args = parser.parse_args()

    arglib.die_if_path_exists( args.output )
    DSSPBINARY = args.dssp_binary
    project = Project.LoadFromHDF( args.project )
    CONF = project.GetEmptyTrajectory()

    run( project, args.output, args.num_procs, 50000 )
