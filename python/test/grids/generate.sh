#!/bin/bash

if [ $# != 1 ]; then
    echo "only got $# arguments"
    echo "Usage: generate.sh [lammps command]";
    exit 0;
fi;

for i in `seq 8`; do
    $1 -var test_index $i < lammps.in;
done;
   

rm COLVAR*
rm log*
