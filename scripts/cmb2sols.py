#!/usr/local/bin/python2.7 -tt

# This script will read a .top file and combine the last two lines if they are both SOL

import re
import sys

inTop = open(sys.argv[1],'r')
outTop = open(sys.argv[1][:-4]+'NEW.top','w')


lines = inTop.readlines()

if lines[-1][0] == '\n':
    lines.pop(-1)

last = lines[-1]
scdLast = lines[-2]

lastSrch = re.search('SOL',last)
scdLastSrch = re.search('SOL',scdLast)

if lastSrch and scdLastSrch:
    lastTmp = last.split()
    scdLastTmp = scdLast.split()

    num = int(lastTmp[1]) + int(scdLastTmp[1])

    newLine = 'SOL             '+str(num)

lines[-2] = newLine
lines.pop(-1)

for line in lines:
    outTop.write(line)


inTop.close()
outTop.close()
