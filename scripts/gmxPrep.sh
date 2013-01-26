#!/bin/bash

# This script will do all the necessary things to prepare
# a $NAME.pdb file for a Gromacs equilibration:
# 1) pdb2gmx
# 2) add_bond
# 3) Edit itp and gro files (modItp.py and modGro.py)
# 4) dstRst 
# 5) BOX (editconf)
BOX="15.77858 6.18062 6.28355"
# 6) emVac
# 7) genbox
# 8) em
# 9) genbox up to 14609
# 10) genion
# 11) em
# 12) make a file containing *.itp, $NAME.top, and the 
#      prepared .gro file for scp'ing to Certainty

function pause() {
    read -p "$*"
}

NAME=$1

echo This script will setup for a GMX equilibration. Using $NAME/$NAME.pdb
echo This script will FAIL IF:
echo "1) amber99.ff is not in the current directory."
echo "2) residuetypes.dat is not in the current directory."
echo "3) modItp.txt must be in the current directory"
echo "4) peptide.gro must be in the current directory"

echo " Do you want to continue? "

pause " Press any key"

cp -r amber99sb.ff $NAME

############################################################
# 1: pdb2gmx
# Make topology files and init.gro


pdb2gmx -f $NAME/$NAME.pdb -o $NAME/init.gro -p $NAME/topol -i $NAME/posre -ff amber99sb -water tip3p -renum -chainsep interactive

echo "------------------------------------------------------"
echo "                add_bond is next"
echo "------------------------------------------------------"
pause "Continue?"

############################################################
# 2: add_bond
# Add tRNA-peptide bond to the topology and remove atoms
# needed to form this bond

add_bond $NAME/topol_Protein.itp $NAME/topol_ProteinNEW.itp $NAME/init.gro $NAME/initNEW.gro 551 713 550 712


# Rename the newly produced file to what it is in topol.top
mv $NAME/topol_Protein.itp $NAME/topol_ProteinOLD.itp
mv $NAME/topol_ProteinNEW.itp $NAME/topol_Protein.itp

# add_bond produces posreNEW.itp in the current dir:
mv $NAME/posre_Protein.itp $NAME/posre_ProteinOLD.itp
mv posreNEW.itp $NAME/posre_Protein.itp
echo "------------------------------------------------------"
echo "                modItp/modGro are next"
echo "------------------------------------------------------"

pause "Continue?"
############################################################
# 3a) modify the itp file to contain the right charges, 
# angles and dihedrals
modItp.py $NAME/topol_Protein.itp modItp.txt

mv $NAME/topol_Protein.itp $NAME/topol_ProteinOLD2.itp
mv $NAME/topol_ProteinNEW.itp $NAME/topol_Protein.itp

# Rename the newly produced .gro file
mv $NAME/init.gro $NAME/initOLD.gro
mv $NAME/initNEW.gro $NAME/init.gro

#############################################################
# 3b) modGro.py
# Change the name of OC2/OC1 -> O in init.gro
modGro.py $NAME/init.gro modGro.txt

echo "------------------------------------------------------"
echo "                   dstRst is next"
echo "------------------------------------------------------"


pause "Continue?"

#############################################################
# 4) Distance Restraints:
# dstRst
# NOTE: It's easier for me to just cd to $NAME and do this...

cd $NAME
echo "In `pwd`"
dstRst.py initNEW.gro 1.5 topol.top ../peptide.gro

cd ..

echo "------------------------------------------------------"
echo "                editconf is next"
echo "NOTE: Choose group 0 (System) when prompted"
echo "------------------------------------------------------"

pause "Continue?"

#############################################################
# 5) BOXING:
# editconf

echo 0 | editconf -f $NAME/initNEW.gro -o $NAME/boxed.gro -c -princ -box $BOX

echo "------------------------------------------------------"
echo "                     emVac is next"
echo "------------------------------------------------------"

pause "Continue?"

#############################################################
# 6) emVac:
# Uses emVac.mdp in the current dir

grompp -f emVac.mdp -c $NAME/boxed.gro -o $NAME/emVac -p $NAME/topol -v

pause "Run this tpr? "

mdrun -s $NAME/emVac -deffnm $NAME/emVac -v


echo "------------------------------------------------------"
echo "                   genbox is next"
echo "------------------------------------------------------"

pause "Continue?"
###########################################################
# 7) Solvation:
# Uses genbox:

genbox -cp $NAME/emVac.gro -cs -p $NAME/topol -o $NAME/slvtd.gro -maxsol 14609

echo "------------------------------------------------------"
echo "                     em is next"
echo "------------------------------------------------------"


pause "Continue?"

############################################################
# 8) em1:
# Uses em.mdp in current dir:

grompp -f em.mdp -c $NAME/slvtd.gro -p $NAME/topol -o $NAME/em -v
pause "Run this tpr? "
mdrun -s $NAME/em -deffnm $NAME/em -v 

echo "------------------------------------------------------"
echo "                   genbox2 is next"
echo "------------------------------------------------------"


##############################################################
# 9) genbox2:
# I need to add up to 14609 total waters

# First I need to figure out how many I added in the first place,
# by reading topol.top

solLine=`tail -n 1 $NAME/topol.top`
solAry=($solLine)
added=${solAry[1]}

toAdd=$(expr 14609 - $added)

if [ $toAdd -ne 0 ]; then
    genbox -cp $NAME/em.gro -cs -p $NAME/topol -o $NAME/slvtd2.gro -maxsol toAdd
else
    cp $NAME/em.gro $NAME/slvtd2.gro
fi

echo "------------------------------------------------------"
echo "                   genion is next"
echo "------------------------------------------------------"


pause "Continue?"

################################################################
# 10) genion
#

grompp -f em.mdp -c $NAME/slvtd2.gro -p $NAME/topol -o $NAME/forGenion

echo 28 | genion -s $NAME/forGenion -p $NAME/topol -o $NAME/neut.gro -nn 1 -np 216

echo "------------------------------------------------------"
echo "                     em2 is next"
echo "------------------------------------------------------"

pause "Continue?"

################################################################
# 11) em

grompp -f em.mdp -c $NAME/neut.gro -p $NAME/topol -o $NAME/em2 -v

pause "Run this tpr file?  "

mdrun -s $NAME/em2 -deffnm $NAME/em2 -v 

################################################################

echo "------------------------------------------------------"
echo "            Placing everything in folder..."
echo "------------------------------------------------------"

mkdir $NAME.toCert

cp $NAME/*.itp $NAME.toCert
cp $NAME/topol.top $NAME.toCert
cp $NAME/em2.gro $NAME.toCert

echo "DONE!"
echo "Folder $NAME.toCert contains everything necessary to continue"
