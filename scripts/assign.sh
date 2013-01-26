#!/bin/bash

echo "Usage: assign.sh Directory FLAGS"
echo "Flags should be contained in quotes"
echo " This script assumes the base directory is one level below the input dir"

cd `dirname $1`

if [[ $1 == "Help" || $1 == "help" || $1 == "--help" || $1 == "-h" ]]
then
  exit
fi


if [ ! -e $1 ]
then
  mkdir $1
fi

if [ ! -e ProjectInfo.h5 ]
then
  echo "No ProjectInfo.h5... Continuing hoping it's in your flags"
fi

if [ ! -e AtomIndices.dat ]
then
  echo "No AtomIndices.dat... Continuing hoping it's in your flags"
fi

Assign.py -g $1/Gens.lh5 -w $1/ $2 &> $1/Assign.log &

wait

