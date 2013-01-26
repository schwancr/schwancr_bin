#!/bin/bash


export GMX_MAXBACKUP=-1
echo "Usage: sh runSims.sh Directory a b" 
echo "Where a and b are the numbers/names of the tprs within directory. "
echo "Note: a.tpr through (b-1).tpr are actually run here"

if [[ $1 == "help" ]]
then
  exit
fi

dir=$1
cd $dir

for (( i=$2; i<$3; i++ ))
do
  while (( `jobs | grep . -c` > 23 ))
  do
    sleep 100
  done
  mdrun -nt 1 -s $i.tpr -cpi $i.cpt -deffnm $i `if [ -e table.xvg ]; then echo "-table table.xvg -tablep table.xvg"; fi` &> /dev/null &
done


wait
