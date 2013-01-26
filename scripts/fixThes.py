#!/usr/bin/env python

import sys

input = open(sys.argv[1],'r')
output = open(sys.argv[1][:-4]+'OUT.pdb','w')


lines = input.readlines()

for line in lines:
	if line[:4] == 'ATOM':
		outLine = line[:70] + '\n'
	else:
		outLine = line

	output.write(outLine)


input.close()
output.close()
