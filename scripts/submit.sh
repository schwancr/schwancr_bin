#!/bin/bash

if [ ! $2 ]
then
  echo "Usage: submit.sh \"Command In Quotes\" JobName queue"
  echo "If queue is blank, the default queue is selected"
  exit
fi

cmd=$1

if [ $3 ]
then
  queue=$3
  if [[ $3 == 'long' || $3 == 'batch' ]]
  then
    walltime='99:00:00:00'
    ppn='24'
  elif [[ $3 == 'SM' ]]
  then
    walltime='99:00:00:00'
    ppn='8'
  elif [[ $3 == 'SP' ]]
  then
    walltime='168:00:00'
    ppn='16'
  elif [[ $3 == 'short_SP' ]]
  then
    walltime='10:00:00'
    ppn='16'
    queue='SP'
  else
    ppn='24'
    walltime='23:59:59'
  fi
else
  queue='default'
  ppn='24'
  walltime='23:59:59'
fi

echo "#!/bin/bash
#PBS -N $2
#PBS -e /home/schwancr/PBS_LOGFILES/$2.\$PBS_JOBID.err
#PBS -o /home/schwancr/PBS_LOGFILES/$2.\$PBS_JOBID.out
#PBS -q $queue
#PBS -l walltime=$walltime
#PBS -l nodes=1:ppn=$ppn

PBS_O_WORKDIR=`pwd`
export PBS_O_WORKDIR

cd \$PBS_O_WORKDIR

$cmd &> LOG.$2.log &

wait

unmountSSHFS.sh" > $2.pbs

qsub $2.pbs

