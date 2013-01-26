#!/bin/bash

# First do vspm24:
fn=/Users/schwancr/DiskUsage.txt

echo "VSPM24:" > $fn
du -h -d 1 /Users/schwancr >> $fn

echo >> $fn
echo "Certainty:" >> $fn
ssh certainty "du -h --max-depth 1 /home/schwancr" >> $fn

echo >> $fn
echo "VSP-COMPUTE-01" >> $fn
ssh vsp-compute-01 "unmountSSHFS.sh; du -h --max-depth 1 /home/schwancr" >> $fn

cat $fn | mail -s "Disk Usages `date`" schwancr@stanford.edu

