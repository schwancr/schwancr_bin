#/bin/bash

if [ ! $# -eq 1 ]
then
  echo "Usage: delete_bad_mem.sh <jobid>"
fi

while ((`checkjob $1 2> /dev/null | grep -c .` ))
do
  node=$( checkjob $1| tr '\n' ' ' | awk -F'Allocated Nodes' '{print $2}' | cut -d[ -f2 | cut -d. -f1 )

  free_MB=`ssh $node "free -m" | grep Mem | awk '{ print $4+$7 }'`

  if [ $free_MB -le 100 ];
  then
    qdel $1
    
    echo "Job $1 killed because memory was $free_MB MB" | mail -s "Job $1" schwancr@stanford.edu
    exit
  fi

  sleep 20

done
