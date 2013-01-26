#!/bin/bash


# This script emulates FAH. It runs frame#.tpr's and then generates the next frame.

# There are four inputs via the command line:
# 1) The directory where the first frame is located
# 2) The first frame#.tpr to start at (i.e. which frame number only)
# 3) The number of ps to extend each frame by
# 4) The number of generations to run

# Note this was originally designed for use with User defined tables, so edit it if you do not
# want to use the -table and -tablep flags in mdrun

echo "Usage: emulateFAH.sh WorkingDirectory FirstFrame ExtendTime(ps) NumberOfGenerations"

cloneDir=$1
tableFile=/home/schwancr/Villin/fahVillin/PROJ0108_files/table.xvg
startGen=$2
extendBy=$3
endGen=$[$2+$4]

cd $cloneDir

for (( i=$startGen; i<$endGen; i++ ))
do

# Do the first run. $CPT is empty so mdrun will run form the stat since it won't find the checkpoint file

   mdrun -s frame$i.tpr -deffnm frame${i} -v -nt 1 -table $tableFile -tablep $tableFile -cpi frame$[i-1].cpt &
   wait

# Now the first frame finished, so generate the next one:

   tpbconv -s frame$i.tpr -extend $extendBy -o frame$[$i+1] &
   wait

done

