#!/bin/bash

if (( ! `ls ~/vspg14c | grep -c .` ))
then
  sshfs server@vspg14c: ~/vspg14c
fi

if (( ! `ls ~/vspm24 | grep -c .` ))
then
  sshfs vspm24: ~/vspm24
fi

if (( ! `ls ~/vspm24_d3 | grep -c .` ))
then
  sshfs vspm24:/Volumes/DiskThree ~/vspm24_d3
fi
