to add this patch to plumed-1.3 it is enough to copy
gromacs-4.6.1 in ./plumed-1.3/patches/diff and 
plumedpatch_gromacs_4.6.1.sh in ./plumed-1.3/patches

Gromacs 4.6 has moved from the former automake configuration
system to CMAKE. To comply with CMAKE we changed the way
in which the patch works.

Example:

tar -zxf gromacs-4.6.1.tar.gz
cd gromacs-4.6.1
export plumedir="PLUMED root"
cp $plumedir/patches/plumedpatch_gromacs_4.6.1.sh .
./plumedpatch_gromacs_4.6.1.sh -patch

now the source code is patched and it is ready
to be 
1. Configured
2. and Compiled

PLUMED works with gromacs configured in any combinations
that doesn't include MULTI-THREADING.
