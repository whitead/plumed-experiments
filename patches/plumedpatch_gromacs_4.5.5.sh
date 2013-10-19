#!/bin/bash
# PATCH SCRIPT FOR GROMACS
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="gromacs-4.5.5"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
RECON_LINKS="$plumedir/recon_src/*.cpp $plumedir/recon_src/*.h"
WHERE_LINKS="./src/kernel/"
# Leave this line here blank for normal patching of plumed
RECON_LIBS= 
# CVS version only : for reconnaissance metadynamics enter here information on the
# locations of the lstdc++ libraries.
#RECON_LIBS="-lstdc++"
# We must define what c++ compiler we are using to compile the reconnaissance routines
RECON_CPP="mpicxx"

PO_FILES="$(
  for file in $plumedir/common_files/*.c
  do f=${file##*/}
  echo $WHERE_LINKS/.deps/${f%.c}.Po
  done
  if [ "$RECON_LIBS" != "" ] ; then 
    for file in $plumedir/recon_src/*.cpp
    do f=${file##*/}
    echo $WHERE_LINKS/.deps/${f%.cpp}.Po
    done
  fi
)"

function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
  touch $PO_FILES
  {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.c}.\$(OBJEXT)"
    done
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo " \\"
        echo -n " ${f%.cpp}.\$(OBJEXT)"
      done
    fi
    echo
    echo -n "PLUMED_SRC="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.c}.c"
    done
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo " \\"
        echo -n " ${f%.cpp}.cpp"
      done
    fi
    echo
    echo
  } > ./src/kernel/plumed.inc
  {
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo "include ./\$(DEPDIR)/${f%.c}.Po"
    done
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo "include ./\$(DEPDIR)/${f%.cpp}.Po"
      done
    fi
  } > ./src/kernel/plumed.Po.inc
  echo "RECON_LIBS=" $RECON_LIBS > ./src/kernel/recon_patch.inc
  if [ "$RECON_LIBS" = "" ] ; then
    { 
      echo "RECON_FLAGS=" 
      echo "RECON_CPP=\$(CC)"
    } >> ./src/kernel/recon_patch.inc
  else
    {
      echo "RECON_FLAGS=-DRECONMETAD" 
      echo "RECON_CPP="$RECON_CPP
    } >> ./src/kernel/recon_patch.inc
  fi
}

function to_do_before_revert () {
  rm $PO_FILES
  rm ./src/kernel/recon_patch.inc
  rm ./src/kernel/plumed.inc
  rm ./src/kernel/plumed.Po.inc
}

function to_do_after_revert () {
  echo > /dev/null
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

