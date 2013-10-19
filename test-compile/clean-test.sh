#!/bin/bash

BASE=`pwd`
cd ..
export plumedir="`pwd`"
cd $BASE
wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-4.6.1.tar.gz
tar -xvf gromacs-4.6.1.tar.gz
cd gromacs-4.6.1
cp ${plumedir}/plumedpatch_gromacs_4.6.1.sh .
./plumedpatch_gromacs_4.6.1.sh -patch
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$BASE/gromacs_test_plumed_mods -DGMX_MPI=ON
make -j 8 && make install
echo "Compiled without error!"
