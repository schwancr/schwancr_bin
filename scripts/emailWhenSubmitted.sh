jobid=$1

state=`checkjob $jobid | grep State | cut -d' ' -f2`

while [[ $state == "Idle" ]]
do
   sleep 60
   state=`checkjob $jobid | grep State | cut -d' ' -f2`
done
 
if [[ $state == "Running" ]]
then
   echo "Job $jobid started at `date`." | mail -s "Job $jobid" schwancr@stanford.edu
else 
   echo "Job $jobid is not Idle, but is also not running..." | mail -s "Job $jobid NOT RUNNING" schwancr@stanford.edu
fi


