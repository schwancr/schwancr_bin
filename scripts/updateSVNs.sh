#!/bin/bash

if [ $# == 0 ]
then
  echo "Usage: updateSVNs.sh PACKAGE_NAME [ PACKAGE_NAME ... ]"
  exit
fi

cd $HOME
export cmd=$(for proj in $*; do echo "cd \$HOME/Installed/$proj; svn update .; python setup.py install"; done)
for machine in vsp-compute-01 certainty
do
  echo $machine
  ssh $machine "source .bash_profile; $cmd; exit"
done

echo "vspm24"
export ARCHFLAGS=""
for proj in $*
do
  cd ~/Installed/$proj
  svn update .
  python setup.py install
done
