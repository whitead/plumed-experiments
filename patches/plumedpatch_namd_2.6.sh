#!/bin/bash
# PATCH SCRIPT FOR NAMD 2.6 
#

# USER DEFINED ARCHITECTURE:
myarch="MacOSX-i686-g++"

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="namd2.6"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
WHERE_LINKS="./src/"

# this changes the names from .c to .C
PLUMED_RENAME_RULE=plumed_rename_rule
function plumed_rename_rule (){
  case "$1" in
  (*.c) echo "${1%.c}.C" ;;
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
  } > ${WHERE_LINKS}/plumed.inc

  cd $myarch ; make depends ; cd ../
}

function to_do_before_revert () {
  rm ${WHERE_LINKS}/plumed.inc
}

function to_do_after_revert () {
  cd $myarch ; make depends ; cd ../
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

