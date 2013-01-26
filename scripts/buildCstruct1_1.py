#!/usr/bin/python
'''

buildCstruct 1.1

Author: Andrea Minoia
Date: March 04th 2010
********************************************
Description:
    Build allotropes structures of carbon, namely graphite sheets and
    carbon nanotubes. Non periodic structure are saturated with hydrogens


Syntax:

    buildCstruct.py [options] outfile

    The available options are:
    --version

        shows program's version number and exit.

    -c, --credits

        display credits.

    -s, --structure

        specify the kind of structure to build. Valid structures are: armcnt, zigzagcnt and hopg.

    -p, --periodicity

        build a periodic structure for TINKER.

    -g, --geometry

        specify the geometry of the structure. For CNTs, use -g index_n cnt_length while
        for hopg use -g size_x size_y (see examples below).
        Size_x, size_y and cnt_length are in Angstom.

    --xyz

        save structure in XYZ format (version 1.1 only).

    --gro

        save structure in gromacs GRO format (version 1.1 only).

    outfile
    is the name of the file where to save the structure.


Requirements:
    Require numpy installed


Known issues and limitations for version 1.1:
    1) Build single wall CNT and rectangular slab of graphtie HOPG only.
    2) Build armchair and zigzag CNT only


Changelog:
    + v 1.1 - March 2010:
        - improved connectivity search algorithm 
        - added support XYZ format
        - added support gromacs GRO format
        
    + v 1.0 - April 2009:
        - first release of buildCstruct


License:
    Freeware. You can use, modify and redistribute the source.
      
    buildCstruct is part of Yasc, a collection of freeware software, therefore it
    cames with _ABSOLUTE_ _NO_(as in NOTHING, NADA, NICTHS, NIENTE, RIEN) WARRANTY.
    The authors of the scripts collected in Yasc are _NOT_ responsible if these software
    will erase your hard disks, empty your bank account, stole your car, seduce
    your wife, shave your dog or make any kind of mess and damages, including loss
    of data or worst. By using Yasc, you _ACCEPT_ these terms.
    Once again: use at your own risk!


Contacts:
    You liked this program? You have suggestions?
    Drop me a mail at minoiaa_at_gmail.com
    Don't forget to visit the wiki: http://chembytes.wikidot.com
'''


#import modules
from sys import argv, exit,path, version

from optparse import OptionParser as OP
try:
    from numpy import zeros, pi, sin, cos, modf, ceil, sqrt
except:
    print "Numpy not installed or not in python path. I give up..."
    exit(10)
from string import lower
from os import path as path_os
from os import system


'''
#===============================================================================
#                               SUBROUTINES
#===============================================================================
'''

def getdist(at1,at2):
    ''' Calculate distance between two particles
    '''
    dist_at=sqrt((at2[0]-at1[0])**2+(at2[1]-at1[1])**2+(at2[2]-at1[2])**2)
    return dist_at

def filecheck(file):
    ''' Check if infile exists
    '''
    if path_os.isfile(file) == 0:
        found=False
    else:
        found=True
    return found

def backup_file(file):
    '''check if file exists. if not open file
    for writing, otherwhise backup the old one in
    #infile_x# with x progressive number
    '''
    tmpvar=file
    count=0
    while 1:
        found=filecheck(file)
        if found:
            count+=1
            file=tmpvar+'_bak-'+str(count)
        else:
            break
    if file !=tmpvar:
        system('mv '+tmpvar+' '+file)

def write_tnk(file,data):
    '''write TINKER xyz/arc
    '''
    file.write(" "+str(len(data))+"\n")
    for line in data:
       outline="%3s  %-3s%12.6f%12.6f%12.6f%6s%6s" % (line[0],line[1],float(line[2])\
                                     ,float(line[3]),float(line[4]),line[5],line[6])
       for i in range(7,len(line)):
          outline=outline+"%6s" % line[i]
       file.write(outline+"\n")
    return

def parsecmd():
    description="Build allotropic structures of Carbon, namely\
graphite hopg and armchair/zigzag carbon nanotubes.\n Output file can\
be saved in TINKER, XYZ or Gromacs GRO formats. Structures can also be periodic.\n"
    usage = "usage: %prog [options] output_file"
    #parse command line
    parser=OP(version='%prog 1.1',description=description, usage=usage)
    
    parser.add_option('-c','--credits',dest='credits',action='store_true',
                     default=False,help='display credits')
    parser.add_option('-s','--struct',dest='structure',default='none',
                      help='define structure: armcnt, zigzagcnt, hopg')
    parser.add_option('-p','--periodic',dest='pbc',action='store_true',
                     default=False, help='build periodic structure for Tinker')
    parser.add_option('-g','--geometry',dest='geometry',nargs=2,type='float',
                      help='define the geometry for the structure')
    parser.add_option('--xyz',dest='xyz',action='store_true',
                      help='write xyz file. This is FAST.')
    parser.add_option('--gro',dest='gro',action='store_true',
                      help='write gromacs gro file. This is FAST too ;).')
    (options, args) = parser.parse_args(argv[1:])
    
    #manage parse errors
    if options.credits: #display credits and quit
        credits="\n**********************************\n\
    Andrea Minoia\n\
    Contacts: minoiaa_at_gmail.com\
              http://chembytes.wikidot.com\
\n*********************************\n"
        print credits
        exit(0)

    if len(args)==0:   #arguments missing
        parser.exit(parser.print_help())
    
    if len(args)>1: #check if more than one argument (NOT OPTION) has been parsed
        parser.error('You have given me more than one argument '+str(args)+'... dunno what to do...\n')
    
    if lower(options.structure) != 'hopg' and lower(options.structure) !='armcnt'\
        and lower(options.structure) !='zigzagcnt':
        parser.error('Uknown structure: valid structures are hopg, armcnt and zigzagcnt')

    return options, args

def armcnt(n,l,ccbond):
    ''' build armchair carbon nanotube
    '''
    atc=[]
    circ1=[]
    circ2=[]
    dx=ccbond*cos(120/2*(pi/180.0))
    dy=ccbond*sin(120/2*(pi/180.0))
    radius=(n*(2*dx+ccbond)+n*ccbond)/(2*pi)
    ycoord=+dy
    natoms=2*n
    #create circumferences
    for i in range(n):
        circ1.append(2*dx+ccbond)
        circ1.append(ccbond)
        circ2.append(ccbond)
        circ2.append(2*dx+ccbond)
    #adjust the circumferences
    circ1.insert(0,0.0)
    circ1.pop()
    circ2.insert(0,dx)
    circ2.pop()
    #Build CNT
    while ycoord>-l:
        ycoord-=dy
        arc=0.0
        for i in range(natoms):
            tmpcoords=['C']
            arc+=circ1[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            atc.append(tmpcoords)
        ycoord-=dy
        arc=0.0
        for i in range(natoms):
            tmpcoords=['C']
            arc+=circ2[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            atc.append(tmpcoords)
    
    pbc_l=abs(ycoord)+dy
    print '\n*******************************'
    print 'armchair CNT: n= ',n,' l (ang)= ',abs(ycoord)
    print 'periodicity (if apply) (ang)= ',pbc_l
    print 'diameter (ang): ',2*radius

    return atc,natoms,pbc_l,len(atc)

def zigzagcnt(n,l,ccbond):
    ''' build armchair carbon nanotube
    '''
    atc=[]
    circ1=[]
    circ2=[]
    dy=ccbond*cos(120/2*(pi/180.0))
    dx=ccbond*sin(120/2*(pi/180.0))
    radius=(n*2*dx)/(2*pi)
    ycoord=+ccbond
    #create circumferences
    for i in range(n):
        circ1.append(2*dx)
    #adjust the circumferences
    
    circ1.pop()
    circ2=list(circ1)  #copy list!!! circ2=circ1 make duplicate and are both modified in the same way
    circ1.insert(0,0.0)
    circ2.insert(0,dx)
    #Build CNT
    while ycoord>-l:
        ycoord-=ccbond
        arc=0.0
        for i in range(n):
            tmpcoords=['C']
            arc+=circ1[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            atc.append(tmpcoords)
        ycoord-=dy
        arc=0.0
        for i in range(n):
            tmpcoords=['C']
            arc+=circ2[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            atc.append(tmpcoords)
        ycoord-=ccbond
        arc=0.0
        for i in range(n):
            tmpcoords=['C']
            arc+=circ2[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            atc.append(tmpcoords)
        ycoord-=dy
        arc=0.0
        for i in range(n):
            tmpcoords=['C']
            arc+=circ1[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            atc.append(tmpcoords)
    pbc_l=abs(ycoord)+ccbond
    print '\n*******************************'
    print 'zigzag CNT: n= ',n,' l (ang)= ',abs(ycoord)
    print 'periodicity (if apply) (ang)= ',pbc_l
    print 'diameter (ang): ',2*radius

    return atc,n,pbc_l,len(atc)

def graphite(x,y,ccbond):
    ''' generate single square sheet of graphite HOPG
    '''
    atc=[]
    dx=ccbond*cos(120/2*(pi/180.0))
    dy=ccbond*sin(120/2*(pi/180.0))           
    ycoord=+dy
    xcoords1=[]
    xcoord=0.00
    xcoords1.append(xcoord)
    #build 1st row for X
    while xcoord<=x:
        xcoord+=ccbond+2*dx
        xcoords1.append(xcoord)
        xcoord+=ccbond
        xcoords1.append(xcoord)
    xcoords1.pop() #remove last element, i.e. the bond exceeding the size
    #build 2nd row for X
    xcoord=dx
    xcoords2=[]
    xcoords2.append(xcoord)
    while xcoord<=x+dx:
        xcoord+=ccbond
        xcoords2.append(xcoord)
        xcoord+=ccbond+2*dx
        xcoords2.append(xcoord)
    xcoords2.pop() #remove last element, i.e. the bond exceeding the size
    while ycoord>-y:
        ycoord-=dy
        for coord in xcoords1:
            tmpcoords=['C']
            tmpcoords.append(coord)
            tmpcoords.append(ycoord)
            tmpcoords.append(0.00)
            atc.append(tmpcoords)
        ycoord-=dy
        for coord in xcoords2:
            tmpcoords=['C']
            tmpcoords.append(coord)
            tmpcoords.append(ycoord)
            tmpcoords.append(0.00)
            atc.append(tmpcoords)
    print '\n*******************************'
    print 'HOPG graphite: a= ',atc[len(xcoords1)-1][1],' b= ',abs(ycoord)
    a_pbc=atc[len(xcoords1)-1][1]+ccbond
    b_pbc=abs(ycoord)+dy
    print 'Periodic (if apply) (ang): a= ',a_pbc, ' b= ',b_pbc

    return atc,len(xcoords1),a_pbc,b_pbc,len(atc)

def connect(coords,natx,pbcx,pbcy,nohcoords):
    '''build connectivity for graphite and nanotube:
    '''
    Ccov_r = 0.77  #covalent radius carbon
    Hcov_r = 0.32  #covalent radius Hydrogen

    btollcc = (2*Ccov_r)*15/100  #bond tollerance of 15%
    btollch = (Hcov_r+Ccov_r)*5/100  #bond tollerance of 5% 
    bondcc=[2*Ccov_r-btollcc,2*Ccov_r+btollcc]
    bondch=[(Hcov_r+Ccov_r)-btollch,(Hcov_r+Ccov_r)+btollch]
    connect=zeros((len(coords),3),int) #init connectivity matrix
    #find connectivity, based on distance
    for i in xrange(len(coords)):
        for j in xrange(i+1,i+2*natx):
            if j<nohcoords:
                at1=[coords[i][1],coords[i][2],coords[i][3]]
                at2=[coords[j][1],coords[j][2],coords[j][3]]
                bond=getdist(at1,at2)
                if bond >= bondcc[0] and bond < bondcc[1]:
                    for k in range(3):
                        if connect[i][k]== 0:
                            connect[i][k]=j+1   #index run from zero, not 1
                            break
                        else:
                            pass
                    for k in range(3):
                        if connect[j][k]== 0:
                            connect[j][k]=i+1   #index run from zero, not 1
                            break
                        else:
                            pass
        if len(coords)!=nohcoords and i<nohcoords: #there are hydrogens
            for j in xrange(len(coords)-nohcoords,len(coords)):
                at1=[coords[i][1],coords[i][2],coords[i][3]]
                at2=[coords[j][1],coords[j][2],coords[j][3]]
                bond=getdist(at1,at2)
                if bond >= bondch[0] and bond < bondch[1]:
                    for k in range(3):
                        if connect[i][k]== 0:
                            connect[i][k]=j+1   #index run from zero, not 1
                            break
                        else:
                            pass
                    for k in range(3):
                        if connect[j][k]== 0:
                            connect[j][k]=i+1   #index run from zero, not 1
                            break
                        else:
                            pass
    #make periodic (TINKER) the structure                
    if pbcx:
        #cicle periodicity along X
        for i in xrange(len(coords)):
            for j in xrange(i,len(coords)):
                if coords[i][1]==0 and j-i==natx-1:
                    for k in range(3):
                        if connect[i][k]== 0:
                            connect[i][k]=j+1   #index run from zero, not 1
                            break
                        else:
                            pass
                    for k in range(3):
                        if connect[j][k]== 0:
                            connect[j][k]=i+1   #index run from zero, not 1
                            break
                        else:
                            pass
    if pbcy:
        i=0
        #cicle periodicity along Y
        while coords[i][2]==0.0:
            j=len(coords)-natx+i
            for k in range(3):
                if connect[i][k]== 0:
                    connect[i][k]=j+1   #index run from zero, not 1
                    break
                else:
                    pass
            for k in range(3):
                if connect[j][k]== 0:
                    connect[j][k]=i+1   #index run from zero, not 1
                    break
                else:
                    pass
            i+=1    
    return connect
            
def data4tnk(coords,conn):
    ''' Prepare all the data to be written in
    tinker format
    '''
    data=[]
    for i in range(len(coords)):
        tmp=[i+1] #progressive index
        for j in range(4): #coordinates
            tmp.append(coords[i][j])
        tmp.append('0') #atom type
        for j in range(3): #connectivity
            if conn[i][j]==0:
                pass
            else:
               tmp.append(conn[i][j])
        data.append(tmp)
    return data

def add_H(coords,natx,structure):
    '''Add hydrogens to nonperiodic structures
    '''
    Hcov_r = 0.32
    Ccov_r = 0.77
    chbond=Hcov_r+Ccov_r
    Hx1=chbond*cos(120/2*pi/180)
    Hy1=chbond*sin(120/2*pi/180)
    if structure =="hopg" or structure=="armcnt": #saturate hopg in y direction or armchair cnt
        for i in xrange(natx):  #upper border
            tmpcoords=['H']
            if modf(float(i)/2.0)[0]==0:
                tmpcoords.append(coords[i][1]+Hx1) #update x even
            else:
                tmpcoords.append(coords[i][1]-Hx1) #update x odd
            if structure=='hopg':   #set y and z for hopg
                tmpcoords.append(coords[i][2]+Hy1)
                tmpcoords.append(0.0)
            else:
                tmpcoords.append(coords[i][2]+Hy1)
                tmpcoords.append(coords[i][3])
            coords.append(tmpcoords)
        for i in xrange(len(coords)-2*natx,len(coords)-natx): #bottom unsaturated border
            tmpcoords=['H']
            if modf(float(i)/2.0)[0]==0:
                tmpcoords.append(coords[i][1]-Hx1)
            else:
                tmpcoords.append(coords[i][1]+Hx1)
            if structure=='hopg':
                tmpcoords.append(coords[i][2]-Hy1)
                tmpcoords.append(0.0)
            else:
                tmpcoords.append(coords[i][2]-Hy1)
                tmpcoords.append(coords[i][3])
            coords.append(tmpcoords)
            
    if structure =="hopg":  #saturate hopg in x direction
        for i in xrange(len(coords)):
            if coords[i][0] == 'C':
                tmpcoords=['H']
                if coords[i][1]==0.0:
                    tmpcoords.append(-chbond)
                    tmpcoords.append(coords[i][2])
                    tmpcoords.append(0.0)
                    coords.append(tmpcoords)
                elif coords[i][1]==coords[natx-1][1]:
                    tmpcoords.append(coords[i][1]+chbond)
                    tmpcoords.append(coords[i][2])
                    tmpcoords.append(0.0)
                    coords.append(tmpcoords)
                    
    if structure == 'zigzagcnt':   #saturate zigzagcnt
        for i in range(len(coords)-natx,len(coords)):
            tmpcoords=['H']
            tmpcoords.append(coords[i][1])
            tmpcoords.append(coords[i][2]-chbond)
            tmpcoords.append(coords[i][3])
            coords.append(tmpcoords)
        for i in range(natx):
            tmpcoords=['H']
            tmpcoords.append(coords[i][1])
            tmpcoords.append(+chbond)
            tmpcoords.append(coords[i][3])
            coords.append(tmpcoords)

def write_xyz(file,data):
    '''
    Write a xyz file.
        Input Variables:
            file: output file (type: file)
            data: list of lists. Each list contains:
                    1) atom name
                2,3,4) X-, Y- and Z-coordinates
        Variables:
            line: store each list of data (type: list)
            outline: string containing a single line to be written in file (type: string)
    '''
    file.write(" "+str(len(data))+"\nGenerated by YASC buildCstruct v1.1\n")
    for line in data:
       outline="%-3s%12.6f%12.6f%12.6f" % (line[0],float(line[1]),float(line[2])\
                                     ,float(line[3]))
       file.write(outline+"\n")
    return

def write_gro(file,data,pbc1="",pbc2=""):
    '''
    Write a gromacs gro file.
        Input Variables:
            file: output file (type: file)
            data: list of lists. Each list contains:
                    1) atom name
                2,3,4) X-, Y- and Z-coordinates
            pbc1/pbc2: periodic lengths 
        Variables:
            line: store each list of data (type: list)
            outline: string containing a single line to be written in file (type: string)
    '''
    file.write("Generated by YASC buildCstruct v1.1\n "+str(len(data))+"\n")
    for index, line in enumerate(data):
       outline="%5i%-5s%5s%5i%8.3f%8.3f%8.3f" % (1,"CNT",line[0],index,float(line[1])/10.0,float(line[2]/10.0)\
                                     ,float(line[3]/10.0))
       file.write(outline+"\n")
    if pbc1 == "":
        outline="  10   10   10\n"
    elif pbc2 == "":
        outline="  10   "+str(float(pbc1)/10.0)+"   10\n"
    else:
        outline="   "+str(pbc1)+"   "+str(pbc2)+"   1\n"
    file.write(outline+"\n")
    return
def main():
    '''
    #===============================================================================
    #                               MAIN MAIN MAIN MAIN
    #===============================================================================
    '''
    ccbond=1.42  #C-C bond length
    
    (options,args)=parsecmd()
    pbcx=False
    pbcy=False

    ofile=args[0] #get output file to save structure
    if lower(options.structure) == "hopg":
        coords,natx,pbc_a,pbc_b,nohcoords=graphite(float(options.geometry[0]),float(options.geometry[1]),ccbond)
        if options.pbc:
            pbcx=True
            pbcy=True
        else:
            add_H(coords,natx,'hopg')
    elif lower(options.structure) == "armcnt":
        coords,natx,pbc_l,nohcoords=armcnt(int(options.geometry[0]),float(options.geometry[1]),ccbond)
        if options.pbc:
            pbcx=False
            pbcy=True
        else:
            add_H(coords,natx,'armcnt')
    else:
        coords,natx,pbc_l,nohcoords=zigzagcnt(int(options.geometry[0]),float(options.geometry[1]),ccbond)
        if options.pbc:
            pbcx=False
            pbcy=True
        else:
            add_H(coords,natx,'zigzagcnt')
    print 'Atoms: ',len(coords)
    print 'saving structure...'
    if not options.xyz and not options.gro:
        conn=connect(coords,natx,pbcx,pbcy,nohcoords) #get connectivity
        tnkdata=data4tnk(coords,conn) #write tinker output file
    backup_file(ofile)
    print '*******************************'
    OUT=open(ofile,'w')
    if options.xyz:
        write_xyz(OUT,coords)
    elif options.gro:
        if options.pbc:
            if lower(options.structure) == "hopg":
                write_gro(OUT,coords,pbc_a, pbc_b)
            else:
                write_gro(OUT,coords,pbc_l)
        else:
                write_gro(OUT,coords)
    else:
        write_tnk(OUT,tnkdata)
    exit(0)
    
if __name__ == "__main__":
    main()
