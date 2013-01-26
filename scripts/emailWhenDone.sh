jobid=$1

state=`checkjob $jobid | egrep "^State:" | cut -d' ' -f2`

while [[ $state == "Running" || $state == "Idle" ]]
do
   sleep 60
   state=`checkjob $jobid | grep "^State:" | cut -d' ' -f2`
done
 
echo "Job $jobid finished at `date`." | mail -s "Job $jobid" schwancr@stanford.edu
