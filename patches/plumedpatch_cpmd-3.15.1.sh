#!/bin/bash
# PATCH SCRIPT FOR CPMD-3.15.1
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# IMPORTANT: modify the lines if necessary!
# ---------------------------------------------------------------------
# a valid C++ compiler:
cxx="mpiCC"
# the corresponding flag to avoid exception handling:
cxxflag="-fno-exceptions"
# on some IBM machines you need one or both the flags in the next line:
# extraflag="-DPLUMED_CPMD_NOUNDERSCORE -DPLUMED_AIXXLF"
extraflag=""
# change this if you want to compile plumed in a different location:
myarch="$PWD"
# ---------------------------------------------------------------------

# definitions specific to this code
CODE="cpmd-3.15.1"
WHERE_LINKS="src-plumed/"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
RECON_LINKS="$plumedir/recon_src/*.cpp $plumedir/recon_src/*.h"
# Leave this line here blank for normal patching of plumed
RECON_LIBS=
# CVS version only : for reconnaissance metadynamics enter here information on the
# location of the lstdc++ libraries (lapack libraries are already linked by CPMD).
#RECON_LIBS="-lstdc++"

# this changes the names from .c to .cpp
PLUMED_RENAME_RULE=plumed_rename_rule
function plumed_rename_rule (){
  case "$1" in
  (*.c) echo "${1%.c}.cpp" ;;
  (*)   echo "$1" ;;
  esac
}

function to_do_before_patch () {
#  echo "Please enter the name of the directory in which you will compile plumed"
#  echo "if it is not the same as the location of the source you are patching i.e. ./"
#  read myarch
#  if [ $myarch == "" ] ; then
#     myarch="."
#  fi
  mkdir src-plumed
}

function to_do_after_patch () {
  {
    echo  "OBJ_PLUMED = \\"
    for file in $plumedir/common_files/*.c
    do
      f=${file##*/}
      echo " ${f%.c}.o \\"
    done
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo " ${f%.cpp}.o \\"
      done
      echo
    fi
    echo  " "
    echo  "SRC_PLUMED = \\"
    for file in $plumedir/common_files/*.c
    do
      f=${file##*/}
      echo " src-plumed/${f%.c}.cpp \\"
    done
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo " src-plumed/${f%.cpp}.o \\"
      done
      echo
    fi
    echo  " "
    echo -n "PLUMED_HEADERS="
    for file in $plumedir/common_files/*.h
      do f=${file##*/}
      echo " \\"
      echo -n " \$(SRC)/${f%.h}.h"
    done
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.h
        do f=${file##*/}
        echo " \\"
        echo -n " \$(SRC)/${f%.h}.h"
      done
    fi
    echo
    echo "RECON_LIBS=" $RECON_LIBS
  } > src-plumed/plumed.inc

  mv $myarch/Makefile $myarch/Makefile.preplumed
  if [ "$RECON_LIBS" != "" ] ; then
    {
      sed "s/CFLAGS =/CFLAGS = -DPLUMED_CPMD -DRECONMETAD ${cxxflag} ${extraflag} /" $myarch/Makefile.preplumed | \
      awk '{print; if($1=="CC") print "CXX =" }' | \
      sed "s/CXX =/CXX = ${cxx}/" | \
      sed 's/OBJECTS =/OBJECTS = $(OBJ_PLUMED) /' | \
      awk '{print; if($1=="OBJ_CC") print "include \$(SRC)/src-plumed/plumed.inc"}' | \
      sed 's/INCFILES =/INCFILES = $(PLUMED_HEADERS) /' | \
      awk '{print; if($1=="\$(CC)"){ 
        printf "\n"
        printf "\$(OBJ_PLUMED) :\n"
        printf "\t\$(CXX) \$(CPPFLAGS) \$(CFLAGS) -c \$(SRC)/src-plumed/\$(@:.o=.cpp)"
      }}' | \
      sed 's/LFLAGS =/LFLAGS = $(RECON_LIBS) /' 
    } > $myarch/Makefile
  else
    {
      sed "s/CFLAGS =/CFLAGS = -DPLUMED_CPMD ${cxxflag} ${extraflag} /" $myarch/Makefile.preplumed | \
      awk '{print; if($1=="CC") print "CXX =" }' | \
      sed "s/CXX =/CXX = ${cxx}/" | \
      sed 's/OBJECTS =/OBJECTS = $(OBJ_PLUMED) /' | \
      awk '{print; if($1=="OBJ_CC") print "include ./src-plumed/plumed.inc"}' | \
      sed 's/INCFILES =/INCFILES = metadyn.h /' | \
      awk '{print; if($1=="\$(CC)"){ 
        printf "\n"
        printf "\$(OBJ_PLUMED) :\n"
        printf "\t\$(CXX) \$(CPPFLAGS) \$(CFLAGS) -c \$(SRC)/src-plumed/\$(@:.o=.cpp)"
      }}'
    } > $myarch/Makefile

  fi

  echo "----------------------------------------------------------------------"
  echo "  NOTE: the Makefile has been modified adding a c++ compiler and"
  echo "  a flag to avoid exception handling."
  echo "  If necessary, please modify variables cxx, cxxflag and extraflag"
  echo "  at the beginning of the patch script plumedpatch_cpmd-3.15.1.sh"
  echo "----------------------------------------------------------------------"
}

function to_do_before_revert () {
  echo > /dev/null  
}

function to_do_after_revert () {
#  echo "Please enter the name of the directory in which you will compile plumed"
#  echo "if it is not the same as the location of the source you are patching i.e. ./"
#  read myarch
#  if [ $myarch == "" ] ; then
#     myarch="."
#  fi
#  mv $myarch/Makefile.preplumed $myarch/Makefile
  rm -fr src-plumed/
}


#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

