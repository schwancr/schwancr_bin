#!/home/schwancr/bin/python -tt


import sys
import re
from numpy import *

''' 
This script reads in a dx file and a pdb, and then determines
the value of the gradient at each point in the pdb.

The strategy is:
1) Find the closest grid point to the point given in the pdb.
2) Determine the slope in the x,y,z direction by determining the
   linear slope in the nr - 1 and nr + 1 points.
3) Output the list of slopes (x,y,z) for the pdb entries
'''
print "Usage: analDx.py pot.dx points.pdb"
print "Output will be potOUT.dx, for whatever filename.dx you gave it"
values = []
n = [ 0,0,0]

def givePt(testPt):
   return values[ n[1] * n[2] * testPt[0] + n[2] * testPt[1] + testPt[2] ]
  
 

def main():
   global n
   global values
   dx = open(sys.argv[1],'r')
   pdb = open(sys.argv[2],'r')
   out = open(sys.argv[1][:-3]+'OUT.txt','w')

   # First, generate the coordinates from the pdb:

   coords = []

   temp = [ 0.,0.,0. ]
   print "Reading pdb coordinates from:",pdb.name
   for line in pdb.readlines():
      if line[:4] == 'ATOM':
         temp[0] = float(line[30:38])
         temp[1] = float(line[38:46])
         temp[2] = float(line[46:54])

         coords.append(list(temp))
       

# Now we have the coordinates to test

#for point in coords:
#   print point

   print "Reading dx values and parameters from:",dx.name

# Now read in the list of values as well as the dx parameters
   values = []
   origin = []
   delta = []
   valueLines = []
   n = [0,0,0]
   dr = ['It',"Didn't",'Work.']
   dxLines = dx.readlines()
   dx.close()
   #count = 0
   N = len(dxLines)
   for line in dxLines:
    #  count+=1
    #  print count, N

      if line[0] == '#':
         continue
      elif line.split()[0] == 'origin':
         origin.append(float(line.split()[1]))
         origin.append(float(line.split()[2]))
         origin.append(float(line.split()[3]))
      elif line.split()[0] == 'delta':
         delta.append(float(line.split()[1])+float(line.split()[2])+float(line.split()[3]))
      elif line.split()[0] == 'object' and line.split()[1] == '1':
         n[0] = int(line.split()[-3])
         n[1] = int(line.split()[-2])
         n[2] = int(line.split()[-1])
      elif len(line.split()) > 3:
         continue
      else:
         valueLines.append(line)

   for line in valueLines:
      temp = line.split()

      for thing in temp:
         values.append(float(thing))
   del dxLines

   dblDelta = [ 2*thing for thing in delta ]

   print "Finding closest grid point and calculating the gradient..."
   count = 0
   for point in coords:
      
      # First find the closest point 
      # Will do this for all three coords
      tempClosePt = []
   
      for i in range(3):
         temp = point[i] - origin[i]
         temp = trunc(temp / delta[i])
         tempClosePt.append( int(temp) )

      #print point, [tempClosePt[i]*delta[i] + origin[i] for i in range(3)] 
      tempClosePt = array(tempClosePt)
      for i in range(3):
         mvArray = array([ (i-2) * (i-1) / 2, i * (i-2) / -1 , i * (i-1) / 2])
         #if (sum(mvArray) != 1):
         #print mvArray
         
         dr[i] = givePt(tempClosePt + mvArray) - givePt(tempClosePt)
         dr[i] = dr[i] / delta[i]
      dr = array(dr)
      dr = dr * 300. * 1.38E-23
      out.write(str(dr[0]) +'\t'+str(dr[1])+'\t'+str(dr[2])+'\n')

      count+=1
   pdb.close()
   out.close()


if __name__ =="__main__":
   main()
