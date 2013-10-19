#!/bin/bash
# PATCH SCRIPT FOR SANDER 11 
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="sander_11"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
WHERE_LINKS="./src/sander/"


function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
#  sed -e "s/^AMBERBUILDFLAGS=/AMBERBUILDFLAGS=-DAMBER /" ${destination}/src/config_amber.h > tmp.h
  echo "AMBERBUILDFLAGS=-DAMBER" >> ${destination}/src/config.h    
  grep -q " -DMPI" ${destination}/src/config.h &&
  {
    sed -e "s/^CPPFLAGS=/CPPFLAGS=-DMPI /" ${destination}/src/config.h > tmp.h
    mv tmp.h ${destination}/src/config.h
  }
  {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n " ${f%.c}.o"
    done
    echo
 } > ${WHERE_LINKS}/plumed.inc
}

function to_do_before_revert () {
  rm ${WHERE_LINKS}/plumed.inc 
}

function to_do_after_revert () {
  sed -e "s/AMBERBUILDFLAGS=-DAMBER//" ${destination}/src/config.h >tmp.h 
  mv tmp.h ${destination}/src/config.h 
  grep -q " -DMPI" ${destination}/src/config.h &&
  {
    sed -e "s/^CPPFLAGS=-DMPI /CPPFLAGS=/" ${destination}/src/config.h > tmp.h
    mv tmp.h ${destination}/src/config.h
  }
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

