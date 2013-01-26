#!/usr/bin/env python

import re
import sys
from subprocess import *


folderList = check_output('ls')


folders = folderList.split('\n')
n = len(folders)
# create list of folders to look in
for (index,folder) in enumerate(reversed(folders) ):
    dotMatch = re.search('\.',folder)
    if dotMatch or folder =='':
        folders.pop(n-index-1)

files = []
# create list of files:
for folder in folders:
    eFile = check_output('ls '+folder+'/full20ns.*.edr',shell=True)
    for file in eFile.split():
        files.append(file)



        
# create .xvg files:
input = open('input.txt','w')

input.write('10')
input.close()
output = open('output.txt','w')

for file in files:
    # First create the output name:
    match = re.search('(\S*)/\S*\.(\d*)\.',file)
    Ename = match.group(1) + '.full20ns.'+match.group(2)
    Popen('g_energy -f '+file+' -o '+Ename+' < '+input.__getattribute__('name'),shell=True,stdout=output)

    #Popen(['g_energy','-f',file,'-o',Ename,'<',input.__getattribute__('name')],stdout=output)
  #  output.seek(0)
  #  numTest = re.search('(\d*)\s*Potential',output.read())
  #  if numTest.group(1) != '12':
  #      print "Potential is not 10!"
  #      exit()
  #  output.seek(0)

output.close()
