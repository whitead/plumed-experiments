#!/bin/bash
# PATCH SCRIPT FOR GROMACS
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

plumedir="/Users/jfdama/Code Repositories/plumed-pinion"

# definitions specific to this code
CODE="gromacs-4.6.1"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
RECON_LINKS="$plumedir/recon_src/*.cpp $plumedir/recon_src/*.h"
WHERE_LINKS="./src/kernel/"
# Leave this line here blank for normal patching of plumed
RECON_LIBS= 
# CVS version only : for reconnaissance metadynamics enter here information on the
# locations of the lstdc++ libraries.
#RECON_LIBS="-lstdc++"
# Leave this line here blank for normal patching of plumed

function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
  {
    echo "set(PLUMED_SRC"
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo "	${f%.c}.c"
    done
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo " ${f%.cpp}.cpp"
      done
    fi
    echo ")"
  } > ./src/kernel/plumed.inc
  if [ "$RECON_LIBS" = "" ] ; then
    { 
      echo "set(RECON_FLAGS )" 
    } > ./src/kernel/recon_patch.inc
  else
    {
      echo "set(RECON_FLAGS -DRECONMETAD)" 
      echo "add_definitions(\${RECON_FLAGS})"
    } > ./src/kernel/recon_patch.inc
  fi
}

function to_do_before_revert () {
  rm ./src/kernel/recon_patch.inc
  rm ./src/kernel/plumed.inc
}

function to_do_after_revert () {
  echo > /dev/null
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

