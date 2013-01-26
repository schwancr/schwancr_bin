#!/home/schwancr/bin/python2.7 -tt

import sys
from numpy import *

finished = False

input = open(sys.argv[1],'r')

name = sys.argv[1][:-4]

lines = input.readlines()

first = lines.pop(0)
second = lines.pop(0)
last = lines.pop(-1)

n = len(lines)

eachN = int( ceil( float(n) / 24.) )

for fileNum in range(24):
    a = fileNum * eachN
    b = (fileNum+1) * eachN
    if b > n:
        finished = True
        b = n
    output = open( name + str(fileNum).zfill(2) + '.gro','w')
    output.write(first)
    output.write(second)

    for line in lines[a:b]:
        output.write(line)

    output.write(last)
    output.close()
    if finished:
        break
input.close()
        
