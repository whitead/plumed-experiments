#!/bin/bash
# PATCH SCRIPT FOR DLPOLY
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="dlpoly_2.20"
WHERE_LINKS="srcmod/Plumed/"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
RECON_LINKS="$plumedir/recon_src/*.cpp $plumedir/recon_src/*.h"
# Leave this line here blank for normal patching of plumed
RECON_LIBS=
# CVS version only : for reconnaissance metadynamics enter here information on the
# locations of the lapack and lstdc++ libraries.
#RECON_LIBS="-lgoto2 -lgfortran -lstdc++"

function to_do_before_patch () {
  mkdir srcmod/Plumed/
}

function to_do_after_patch () {
  # to be removed whem a pbc routine common to all the codes will be written
  awk '/subroutine images/{print;for(i=0;i<230;i++){getline;print;if($2=="subroutine")exit}}' \
       srcmod/utility_module.f >  srcmod/images.f
  awk '/subroutine invert/{print;for(i=0;i<230;i++){getline;print;if($2=="subroutine")exit}}' \
       srcmod/utility_module.f >> srcmod/images.f
 {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	Plumed/${f%.c}.o"
    done
    echo 
    echo -n "OBJ_RECON="
    if [ "$RECON_LIBS" != "" ] ; then 
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo " \\"
        echo -n " Plumed/${f%.cpp}.o"
      done
    fi
    echo
    echo -n "HEAD_RECON="
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.h
        do f=${file##*/}
        echo " \\"
        echo -n " Plumed/${f%.h}.h"
      done
    fi
    echo
    echo "RECON_LIBS="$RECON_LIBS
    if [ "$RECON_LIBS" = "" ] ; then
      echo -n "RECON_FLAGS=" 
    else
      echo -n "RECON_FLAGS=-DRECONMETAD" 
    fi
  } > ./srcmod/Plumed/plumed.inc
}

function to_do_before_revert () {
  rm srcmod/images.f
}

function to_do_after_revert () {
  rm -fr srcmod/Plumed/
}


#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

