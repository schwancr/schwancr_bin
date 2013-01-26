#!/bin/bash
jobid=$1
pbsFN=$2

state=`checkjob $jobid | grep State | cut -d' ' -f2`

while [[ $state == "Idle" ]]
do
   sleep 60
   state=`checkjob $jobid | grep State | cut -d' ' -f2`
done

while [[ $state == "Running" ]]
do
   sleep 60
   state=`checkjob $jobid | grep State | cut -d' ' -f2`
done

newJobID=`qsub $2`

emailWhenSubmitted.sh $newJobID
