#!/bin/bash

# this script will get a pdb from the database via ftp, then analyze it with seqCO.py
# then it will delete the pdb file.


cat $1 | while read line 
do
   export file=${line}

# first get the file:

   export midTwo=`echo $file | cut -c2,3 | tr "[:upper:]" "[:lower:]"`

   export name=`echo $file | tr "[:upper:]" "[:lower:]"`


   ftp -V ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/$midTwo/pdb$name.ent.gz . 

# unzip the file
   gzip -d pdb$name.ent.gz

# rename and run seqCO
   mv pdb$name.ent $name.pdb

   seqCO.py $name.pdb

   wait

# remove the pdb once the calculation is done
   rm $name.pdb

done
