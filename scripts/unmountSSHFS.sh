#!/bin/bash

#DIRS="/home/schwancr/vspg14c /home/schwancr/vspm24 /home/schwancr/vspm24_d3 /home/schwancr/vspg13b /home/schwancr/cert"
#DIRS="/home/schwancr/vspm24 /home/schwancr/cert /home/schwancr/cert-shared /home/schwancr/vspm24_vol /home/schwancr/vspg14c /home/schwancr/vspg13b"
#DIRS="/home/schwancr/vspm24 /home/schwancr/vspm24_vol /home/schwancr/vspg14b /home/schwancr/vspm44"
DIRS="/home/schwancr/vspm24 /home/schwancr/vspm24_vol"

for dir in $DIRS;
do
   #if (( `ls $dir | grep -c .` ))
   #then
   #   fusermount -u -z $dir
   #fi
   fusermount -uz $dir
done
