#!/bin/bash

if [ -e $1 ]
then
	echo "$1 already exists!"
	exit
fi

echo "#!/usr/bin/env python" >> $1
echo " " >> $1
echo "from argparse import ArgumentParser" >> $1
echo "parser = ArgumentParser()" >> $1
echo "parser.add_argument('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')" >> $1
echo "parser.add_argument('-a',dest='ass_FN',default='./Assignments.Fixed.h5',help='Assignments from msmbuilder [ ./Assignments.Fixed.h5 ]')" >> $1
echo " " >> $1
echo "args = parser.parse_args()" >> $1
echo " " >> $1
echo "import numpy as np" >> $1
echo "from msmbuilder import Serializer" >> $1
echo "import matplotlib" >> $1
echo "matplotlib.use('pdf')" >> $1
echo "from matplotlib.pyplot import *" >> $1
echo "from pyschwancr import dataIO" >> $1
echo "import os, sys, re" >> $1
echo " " >> $1

vi $1
