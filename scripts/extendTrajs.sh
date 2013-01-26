

if [ -e $1/ALREADY_EXTENDED_1 ]
then
  echo "Already done. exiting"
  exit
fi

export GMX_MAXBACKUP=-1

numGRO=`ls $1 | grep gro -c`
num=`ls $1 | grep tpr -c`

if [[ $numGRO < $num ]]
then
  echo "Not all gro's are present. Exiting"
  exit
fi


for (( i=0; i<$num; i ++ ))
do
  mv $1/$i.tpr $1/${i}_old.tpr
  tpbconv -s $1/${i}_old.tpr -extend 1000 -o $1/$i.tpr
done


echo "Extended this directory" > $1/ALREADY_EXTENDED_1
