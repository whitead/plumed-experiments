#!/bin/bash
# PATCH SCRIPT FOR NAMD 2.7 
#

# USER DEFINED ARCHITECTURE:
myarch="Linux-x86_64-g++"

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="namd2.7"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
RECON_LINKS="$plumedir/recon_src/*.cpp $plumedir/recon_src/*.h"
WHERE_LINKS="./src/"
# Leave this line here blank for normal patching of plumed
RECON_LIBS=
# CVS version only : for reconnaissance metadynamics enter here information on the
# locations of the lapack and lstdc++ libraries.
# RECON_LIBS="-llapack -lstdc++"

# this changes the names from .c to .C
PLUMED_RENAME_RULE=plumed_rename_rule
function plumed_rename_rule (){
  case "$1" in
  (*.c) echo "${1%.c}.C" ;;
  (*)   echo "$1" ;;
  esac
}

# this changes the names from .cpp to .C
RECON_RENAME_RULE=recon_rename_rule
function recon_rename_rule (){
  case "$1" in
  (*.cpp) echo "${1%.cpp}.C" ;;
  (*)   echo "$1" ;;
  esac
}

function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
 {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo "	\\"
      echo -n "	\$(DSTDIR)/${f%.c}.o"
    done
    echo
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo "  \\"
        echo -n " \$(DSTDIR)/${f%.cpp}.o"
      done
      echo " "
      echo "RECON_FLAGS=\$(COPTD)RECONMETAD" 
    else
      echo " "
      echo "RECON_FLAGS=" 
    fi
    echo "RECON_LIBS="$RECON_LIBS
  } > $WHERE_LINKS/plumed.inc

  cd $myarch ; ln -s $WHERE_LINKS/plumed.inc .; make depends ; cd ../
}

function to_do_before_revert () {
  rm $myarch/plumed.inc; rm $WHERE_LINKS/plumed.inc
}

function to_do_after_revert () {
  cd $myarch ; make depends ; cd ../
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

