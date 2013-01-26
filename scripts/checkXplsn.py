#!/usr/bin/env python

import sys
import numpy

'''
This script will determine whether any amino acids
have exploded. It will do this by testing the distance
between atoms in all residues.

The script reads in a gro file (ex: input.gro) and then outputs to 
input.errors.txt the residues that did not satisfy a given
criterion
'''
critDist = 1.5 # nm

skipList = ['SOL','HOH','CL ','NA ','MG ','CD ','K  ']

outTxt = open(sys.argv[1][:-3]+'errors.txt','w')
inGro = open(sys.argv[1],'r')
def main():



    print sys.argv[0]

    lines = inGro.readlines()

    lines.pop(-1)
    lines.pop(0)
    lines.pop(0)

    lastLine = 11
    resLines = []

    atms = len(lines)

    for index,line in enumerate(lines):
#        print int(line[0:5]), lastLine
        print numpy.ceil( float(index) * 100. / float(atms) ), '%\r',
        if line[5:8] in skipList:
            lastLine = int(line[0:5])
            continue

        if int(line[0:5]) == lastLine:
            resLines.append(line)
        else:
            maxDist(resLines)
            resLines = []
            resLines.append(line)


        lastLine = int(line[0:5])
#        print 'lastLine: ',lastLine,
#        print 'int(line[0:5]: ',int(line[0:5])

    inGro.close()
    outTxt.close()

    print '\n   Completed. Output printed to ', sys.argv[1][:-3]+'errors.txt\n'

def maxDist(atomLines):
    
    coords = []

    for line in atomLines:
        x = float(line[20:29])
        y = float(line[29:36])
        z = float(line[37:44])

        coords.append( numpy.array((x,y,z)) )

    testDist = 0.

    n = len(coords)

    failed = False

    for i in range(n):
        if failed == True:
            break
        for j in range(n):
            if j==i:
                continue

            diff = coords[i] - coords[j]
            testDist = numpy.sqrt( numpy.dot( diff, diff ) )
            
            if testDist > critDist:
                failed = True
                outStr = 'Residue '+atomLines[0][:8]+' exploded! (atoms: '+atomLines[0][15:20]+' to '+atomLines[n-1][15:20]+')\n'
                
                outTxt.write(outStr)
                print outStr
                break


    return

if __name__=='__main__':
    main()
