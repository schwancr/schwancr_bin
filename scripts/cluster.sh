#!/bin/bash

echo "Usage: cluster.sh Directory FLAGS"
echo "Flags should be contained in quotes"
echo "   Here are your options:"
echo "      -k # of state"
echo "      -u Stride to sample for clustering"
echo "      -r Radius for states"
echo "      -m Local K-Medoid iterations"
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

Cluster.py -w $1/ $2 &> $1/Cluster.log &

wait

