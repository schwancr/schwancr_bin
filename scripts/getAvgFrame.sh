#!/bin/bash

# This script outputs the run# followed by the average frame in the given clones.

runList=`ls | grep RUN`

for run in $runList
do
  cloneList=`ls $run | grep CLONE`
  sum=0
  count=0
  for clone in $cloneList
  do
    count=$[ $count + 1 ]
    sum=$[ $sum + `ls $run/$clone | grep -c xtc` ] 
  done
  echo ${run}:  $[ $sum / $count ] xtc\'s present on average in clones.
done
