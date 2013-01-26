#!/bin/bash


qnodes -x | sed -e "s/<Node>/\n\n<Node>/g" | grep "<state>free</state>" | cut -d">" -f3 | cut -d"<" -f1


