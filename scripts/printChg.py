#!/usr/local/bin/python2.7 -tt

import re
import sys

# This script inputs a topology file and a gro file
# and outputs a modified gro file with an extra entry
# corresponding to the partial charge located in the
# topology.


if len(sys.argv) < 3 or sys.argv[1]=='help':
    print "\n\nThis function will input two files, and output one file"
    print "\n  Input 1: conf.gro; This file is the gro file with the"
    print "\tatom coordinates."
    print "\n  Input 2: topol.top; This file is the topology file with"
    print "\tpartial charges associated with the atoms in conf.gro"
    print "\tNOTE: This is the topol.top which points to topol_*.itp files"
    print "\n  Optional Input 3: Name of the confMod.gro file with partial"
    print "\tcharges included. If blank the default is input2Mod.gro"
    print "\n  Output is the same gro file but with an added column that has"
    print "\tthe atom's partial charge\n"
    exit()

if len(sys.argv) != 4:
    outGroM = open(sys.argv[1][:-4]+'Mod.gro','w')
else:
    outGroM = open(sys.argv[3],'w')
inGro = open(sys.argv[1],'r')
inTop = open(sys.argv[2],'r')


# First we must generate the list of partial charges
# by iterating through the .itp files

charges = []

itpFiles = []
print "Reading topology ("+inTop.name+")..."
for topLine in inTop:
    matchObj = re.search(r'"(\w+\.itp)"',topLine)
    if matchObj:
        itpFiles.append( matchObj.group(1) )

for itp in itpFiles:
    itpFile = open(itp,'r')
    atmLast = 0
    for itpLine in itpFile:
        matchObj = re.search('\; qtot|\;  CRS',itpLine)
        if matchObj:
            # Then we have an atom line.
            charges.append( itpLine.split()[6] )
            atmNew = int(itpLine.split()[0])
            if atmNew != atmLast + 1:
                print 'Atom missed before this one!: ',itpLine
            atmLast = atmNew

groLines = inGro.readlines()

#print len(charges)

XSlines = ['','','']

XSlines[2] = groLines.pop(-1)
XSlines[0] = groLines.pop(0)
XSlines[1] = groLines.pop(0)

n = len(charges)

otherCharges = { 'CL':'-1.00',
                 'MG':'2.000',
                 'NA':'1.000',
                 'K' :'1.000',
                 'HW1':'0.417',
                 'HW2':'0.417',
                 'OW':'-0.834' }

# The HW and OW charges are based on TIP3P in amber99.ff dist
# in GROMACS
print "Writing modified gro file ("+outGroM.name+")..."

outGroM.write(XSlines[0])
outGroM.write(XSlines[1])

for index,line in enumerate(groLines):
    if index < n:
        # Then we have the atom's charge stored in charges
        atmCharge = ' ' + charges[index] + '\n'
    else:           # Then we need to generate the charge
        # This should only be one of: CL, MG, NA, K, HW, OW
        atmName = line[8:15].strip()
        if atmName not in otherCharges:
            print 'Atom '+atmName+' is not in the dictionary defined in this script!'
            exit()
        atmCharge = ' ' + otherCharges[ atmName ] + '\n'
        
    
    printLine = line[:-1] + atmCharge
    
    outGroM.write(printLine)

outGroM.write(XSlines[2])

inGro.close()
inTop.close()
outGroM.close()

print "Done! \nOutput saved to: "+outGroM.name+'\n'
