#!/bin/bash

a=( $( echo $( ls $1/Sims | grep trj ) ) )

numDirs=${#a[*]}
crtDir=`pwd`
if (( $numDirs <= 50 ))
then
  max=1
else
  max=$numDirs
fi

for (( i=0; i < $max; i+=50 ))
do
  trjList=$(
  for (( j=$i; j< $i+50; j++ ))
  do
     echo $1/Sims/${a[$j]}
  done)

  submit.sh "for i in $( echo $trjList ); do runSims.sh \$i 0 100; done" ${1}_pfoldSims_$i `if [ $max != 1 ]; then echo long; fi`
done
