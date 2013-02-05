#!/bin/bash

#MACHINES=(vspg13b vspg14c vspm24 vspm24 certainty)
#DIRS=("" "" "" /Volumes/DiskThree) # "" Means use the default which is the remote home directory.
#TARGETS=(~/vspg13b ~/vspg14c ~/vspm24 ~/vspm24_d3 ~/cert)
#USERS=(server server schwancr schwancr schwancr)

MACHINES=(vspm24 vspm24)
DIRS=("" /Volumes)
TARGETS=(~/vspm24 ~/vspm24_vol)
USERS=(schwancr schwancr)

#MACHINES=(vspm24 vspm24 vsp12a vsp12b)
#DIRS=("" /Volumes "server2/data/SVR1716710858" "server2/data/SVR1716710859")
#TARGETS=(~/vspm24 ~/vspm24_vol ~/vsp12a ~/vsp12b)
#USERS=(schwancr schwancr server server)

#MACHINES=(vspm24 vspm24 vspg14b vspm44)
#DIRS=("" /Volumes "server2/data/SVR1716465104" /Volumes)
#TARGETS=(~/vspm24 ~/vspm24_vol ~/vspg14b ~/vspm44)
#USERS=(schwancr schwancr server kai)

#MACHINES=(vspm24 vspm24 certainty certainty vspg14c)
#DIRS=("" /Volumes "/panfs/panfs-certainty/home-pande/shared" "" "")
#TARGETS=(~/vspm24 ~/vspm24_vol ~/cert-shared ~/cert ~/vspg14c)
#USERS=(schwancr schwancr schwancr schwancr server)

#if (( ! `screen -list | grep -c vspcert` ))
#then
#  sh ~/schwancr_bin/certLogin.sh
#fi

for (( i=0; i < ${#MACHINES[*]}; i++ ))
do
  host=${MACHINES[$i]}
  target=${TARGETS[$i]}
  dir=${DIRS[$i]}
  user=${USERS[$i]}
  #echo $host 
  if [ ! -e $target ]
  then
    mkdir $target
  fi
  
  if (( ! `ls $target | grep -c .` ))
  then
    sshfs ${user}@${host}:$dir $target -o transform_symlinks
  fi
done

