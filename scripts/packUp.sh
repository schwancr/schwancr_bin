#!/bin/bash



cd $1.toCert

mkdir $1.packUp

pack=$1.packUp

cp -i topol.top $pack
cp -i *.itp $pack
cp -ir $1 $pack
cp -i full*.gro $pack


if [ ! -e ../packed ]; then
  mkdir ../packed
fi

mv $pack ../packed

cd ..
