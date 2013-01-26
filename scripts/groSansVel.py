#!/home/schwancr/bin/python2.7 -tt

import sys
import re

def main():

    groIn = open(sys.argv[1],'r')
    groOut = open(sys.argv[1][:-4]+'_sansVel'+sys.argv[1][-4:],'w')

    line = groIn.readline()
    groOut.write(line)
    line = groIn.readline()
    groOut.write(line)
    for line in groIn:
        groOut.write(line[:44]+'\n')

    groOut.close()
    groIn.close()
if __name__ == '__main__':
    main()
