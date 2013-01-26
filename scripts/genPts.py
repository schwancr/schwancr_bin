#!/usr/bin/python2.6 -tt

import sys
from numpy import *
import re
from random import *
'''
This script will generate a grid of points given two initial points.
It will draw a line between the given points, then produce a vector 
perpindicular to that, and then another vector perpindicular to both.

This defines the grid system, and points are generated with an inputted
spacing and size.
'''

print "Usage: genPts.py input.txt"


if sys.argv[1].lower() in ['help','--help','-h','-help','--h']:
   print "This script will generate a grid of points given two initial points."
   print "It will draw a line between the given points, then produce a vector"
   print "perpindicular to that, and then another vector perpindicular to both."

   print "\nThis defines the grid system, and points are generated with an inputted"
   print "spacing and size."

   print "The input.txt file has five variables to set:"
   print "1) r0"
   print "2) r1"
   print "3) Nray (number of points to add perpindicular to the line"
   print "4) Spacing for the ray points"
   print "5) Nline (number of points to generate along the line"
   
   exit()


##### BEGIN:

# First read in all the data:

input = open(sys.argv[1],'r')

temp = input.readline().split()

r0 = []
r1 = []

for thing in temp:
   r0.append(float(thing))
   
r0 = array(r0)

temp = input.readline().split()

for thing in temp:
   r1.append(float(thing))
   
r1 = array(r1)

Nray = int(input.readline())
spacing = float(input.readline())
Nline = int(input.readline())

# Next generate the vector corresponding to the line i

v = r1 - r0

w = array([uniform(-1,1),uniform(-1,1),0.])

w[2] = - ( v[0]*w[0] + v[1]*w[1] ) / v[2]

# Now I have v,w that define my coordinate system.

x = array([ uniform(-1,1),0.,0. ])

A = array([[v[1],v[2]],[w[1],w[2]]])

C = array([-v[0]*x[0],-w[0]*x[0]])

x[1],x[2] = dot(linalg.inv(A),C)

x = x / sqrt(dot(x,x)) * spacing
w = w / sqrt(dot(w,w)) * spacing
v = v / sqrt(dot(v,v)) * sqrt(dot(r1-r0,r1-r0)) / float(Nline)

# Now we have the coordinates to make the grid points. 

start = r0 - Nray * x - Nray * w

out = open(sys.argv[1][:-4]+'OUT.pdb','w')
count=0
for ix in range(Nray * 2 + 1):
   for iw in range(Nray * 2 + 1):
      for iv in range(Nline + 1):
         count+=1
         temp = start + x * ix + w * iw + v * iv
         outLine = "ATOM  "+str(count).rjust(5,)+"  CA"+' '+' ALA'+' '+'A'+'   1'+' '+'   '+"%(x)8.3f%(y)8.3f%(z)8.3f" % {'x':temp[0],'y':temp[1],'z':temp[2]} + "  1.00  0.00           C \n"
         #out.write(str(temp[0])+'\t'+str(temp[1])+'\t' +str(temp[2])+'\n')
         out.write(outLine)
print "This generated ", (2*Nray+1)**2*(Nline+1), " points."
out.close()
input.close()
