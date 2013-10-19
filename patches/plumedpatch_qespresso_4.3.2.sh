#!/bin/bash
# PATCH SCRIPT FOR GROMACS
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="qespresso-4.3.2"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
RECON_LINKS="$plumedir/recon_src/*.cpp $plumedir/recon_src/*.h"
WHERE_LINKS="./clib/"
# Leave this line here blank for normal patching of plumed
RECON_LIBS=
# CVS version only : for reconnaissance metadynamics enter here information on the
# locations of the lapack and lstdc++ libraries.
#RECON_LIBS="-lgoto2 -lgfortran -lstdc++"
# CVS version only : Here you must define a c++ compiler that can be used to compile 
# the recon metadynamics files 
#RECON_CPP="g++"


function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
  {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.c}.o"
    done
if [ "$RECON_LIBS" != "" ] ; then
    for file in $plumedir/recon_src/*.cpp
      do f=${file##*/}
      echo " \\"
      echo -n " ${f%.cpp}.o"
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
  } > ./clib/plumed.inc
  mv make.sys make.sys.original
  if [ "$RECON_LIBS" = "" ] ; then 
    awk '{if($1=="DFLAGS")print $0,"-DPLUMED_QESPRESSO"; else print}' make.sys.original > make.sys
  else
    awk -v recon="$RECON_LIBS" -v comp="$RECON_CPP" '{if($1=="DFLAGS")print $0,"-DPLUMED_QESPRESSO -DRECONMETAD"; else if( $1=="\$(CC)" ) print $0 "\n \n%.o: %.cpp \n\t"comp," $(CFLAGS)  -c $<"; else if( $1=="LIBS" ) print "RECON_LIBS  =", recon,"\n",$0,"$(RECON_LIBS)"; else print}' make.sys.original > make.sys
  fi
}

function to_do_before_revert () {
  rm ./clib/plumed.inc
}

function to_do_after_revert () {
  echo > /dev/null
  mv make.sys.original make.sys
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

