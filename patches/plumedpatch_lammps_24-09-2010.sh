#!/bin/bash
# PATCH SCRIPT FOR LAMMPS 
#
# this script needs to be launched from the src directory in lammps 
#
#
destination="$PWD/src"

# definitions specific to this code
CODE="lammps_24-09-2010"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
WHERE_LINKS="./USER-PLUMED"

# this changes the names from .c to .C
PLUMED_RENAME_RULE=plumed_rename_rule
function plumed_rename_rule (){
  case "$1" in
  (*.c) echo "${1%.c}.cpp" ;;
  (*)   echo "$1" ;;
  esac
}

function to_do_before_patch () {
  echo > /dev/null
  mkdir ${destination}/USER-PLUMED
}

function to_do_after_patch () {
#
# after creating links and patches
#
cd USER-PLUMED
cat > ./Install.sh << \EOF
# Install/unInstall package classes in LAMMPS

if (test $1 = 1) then

  for file in *.cpp *h;do
    cp -p $file ../
  done

elif (test $1 = 0) then

  for file in *.cpp *h;do
    rm ../$file
  done

fi
EOF
cat > ./Install.csh << \EOF
# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  foreach file (*cpp *.h)
    cp $file ../
  end

else if ($1 == 0) then

  foreach file (*cpp *.h)
    rm ../$file
  end

endif
EOF
cat >./fix_plumed.cpp << \EOF
/*
*******************************************************************************
*                                                                             *
*                                PLUMED                                       *
*   A Portable Plugin for Free Energy Calculations with Molecular Dynamics    *
*                              VERSION 1.2                                    *
*                                                                             *
*******************************************************************************
*
*  
*  Copyright (c) 2010 The PLUMED team.
*  See http://www.plumed-code.org for more information. 
*
*  This file is part of PLUMED.
*
*  PLUMED is free software: you can redistribute it and/or modify
*  it under the terms of the GNU Lesser General Public License as 
*  published by the Free Software Foundation, either version 3 of 
*  the License, or (at your option) any later version.
*
*  PLUMED is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General
*  Public License along with PLUMED.  
*  If not, see <http://www.gnu.org/licenses/>.
*
*  For more info, see:  http://www.plumed-code.org
*  or subscribe to plumed-users@googlegroups.com
*
*/


#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "fix_plumed.h"

static const double PLUMED_UNSET=-1.23456789e12;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixPlumed::FixPlumed(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  double *mss,*chg;
  int *int_p,size;
  int i,j,k,i_c,nn,mm;
  int me;

  if (narg != 7) error->all("Illegal fix plumed command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&size);

  if (!atom->tag_enable)
    error->all("fix plumed requires atom tags");

  if (!atom->tag_consecutive())
    error->all("fix plumed requires consecutive atom tags");

  if (strcmp(arg[1],"all") != 0)
    error->all("fix plumed requires using group all");

//  printf("PROCESSOR %d SIZE %d\n",me,size);

  if (strcmp(arg[3],"plumedfile") == 0) {
     printf("PLUMED: FOUND THE INPUT FILE %s\n",arg[4]);
     if (strcmp(arg[5],"outfile") == 0) {
     printf("PLUMED: FOUND THE OUTPUT FILE %s\n",arg[6]);
     } else {
          error->all("Illegal fix plumed command");
     }
     // detect the number of atoms, mass and charge
     // number of atoms in this group (default must be all )
     nn=int(atom->natoms);
     atom->check_mass();
  printf("PROCESSOR %d NATOMS %d\n",me,nn);
     mss=(double *)malloc(nn*sizeof(double));          
     chg=(double *)malloc(nn*sizeof(double));          
  //   printf("FOUND NATOMS FOR GROUP %d : %d\n",0,nn);
  //   printf("FOUND MASS  %f DT %f\n",group->mass(0),update->dt);
     //for(i=0;i<nn;i++){chg[i]=atom->q[i];} 
     for(i=0;i<nn;i++){chg[i]=PLUMED_UNSET;} 
     for(i=0;i<nn;i++){mss[i]=PLUMED_UNSET;} 
     //use the groups for massses
     MPI_Request request;
     int groupbit = group->bitmask[igroup]; //get the groupbit of all atoms
     // retrieve masses 
     if (atom->rmass) {
      for (int i = 0; i < atom->nlocal; i++){
        if (atom->mask[i] & groupbit) { 
             mm=atom->tag[i]-1; 
             mss[mm] = atom->rmass[i];
//             printf("I AM PROC %d AND I OWN AT %d R-MASS %f TAG %d \n",me,mm,mss[mm],atom->tag[i]);
             if (me) MPI_Isend(&(mss[mm]),1,MPI_DOUBLE,0,mm,world,&request);
        }
      }
     } else {
        for (int i = 0; i < atom->nlocal; i++){
          if (atom->mask[i] & groupbit){ 
               mm=atom->tag[i]-1; 
               mss[mm] = atom->mass[atom->type[i]];
//               printf("I AM PROC %d AND I OWN AT %d MASS %f TAG %d \n",me,mm,mss[mm],atom->tag[i]);
               if (me) MPI_Isend(&(mss[mm]),1,MPI_DOUBLE,0,mm,world,&request);
          }
        }
     }

//     printf("ME: %dEVERYONE SENT \n",me); 

     if(!me){
        MPI_Status status;
        for (i=0;i<nn;i++){
          if(mss[i]==PLUMED_UNSET)MPI_Recv(&(mss[i]),1,MPI_DOUBLE,MPI_ANY_SOURCE,i,world,&status); 
        }
     }

//     printf("ME: %d DONE RECEIVE : Q_FLAG %d \n",me,atom->q_flag); 

     // retrieve charges

     for (i = 0; i < atom->nlocal; i++){
        if (atom->mask[i] & groupbit){ 
             mm=atom->tag[i]-1; 
            if(atom->q_flag){
               chg[mm] = atom->q[i];
             } else  {
               chg[mm] = 0.0 ;
             }
//             printf("I AM PROC %d AND I OWN AT %d CHG %f \n",me,mm,chg[mm]);
             if (me) MPI_Isend(&(chg[mm]),1,MPI_DOUBLE,0,mm,world,&request);
        }
     }

//     printf("ME: %dEVERYONE SENT \n",me); 

     if(!me){
        MPI_Status status;
        for (i=0;i<nn;i++){
          if(chg[i]==PLUMED_UNSET)MPI_Recv(&(chg[i]),1,MPI_DOUBLE,MPI_ANY_SOURCE,i,world,&status); 
        }
     }
   
//     printf("ME: %d DONE RECEIVE  \n",me); 

     MPI_Bcast(mss,nn,MPI_DOUBLE,0,world);
     MPI_Bcast(chg,nn,MPI_DOUBLE,0,world);
  
//     if(!me){for(i=0;i<nn;i++)printf("ME %d I %d MASS %f CHG %f\n",me,i,mss[i],chg[i]);}
      //  the plumed object 
     int_p=(int *)malloc(atom->natoms*sizeof(int));
     for(i=0;i<atom->natoms;i++)int_p[i]=0;
      
     if(!me){
         plumed_ptr= new Plumed(arg[4],arg[6],&nn,mss,chg,&(update->dt),force->boltz); 
         // delete the vectors
         // create a group with only the needed atoms
         int narg_p;                 
         // now create a  group for plumed purpose 
         for(i_c=0;i_c<plumed_ptr->colvar.nconst;i_c++){
              for(j=0;j<plumed_ptr->colvar.natoms[i_c];j++){
                  int_p[plumed_ptr->colvar.cvatoms[i_c][j]]=1; // int p has the "c" like definition (start from 0 )    
                  //printf("PROC %d ATOMS: CV %d PROG %d IND %d \n",me,i_c,j,plumed_ptr->colvar.cvatoms[i_c][j]);
              }
         }
     }
     MPI_Bcast(int_p,nn,MPI_INT,0,world);

     group->create_plumed("plumed",int_p);
     plumedgroup=group->find("plumed");
     if (plumedgroup == -1) 
       error->all("Fix plumed group ID does not exist"); 
      plumedgroup2bit = group->bitmask[plumedgroup];
     //for (int iarg = 0; iarg < narg_p; iarg++) { fprintf(screen,"ARG %s ",arg_p[iarg]);}fprintf(screen,"\n");
     //printf("ASSIGNMENT DONE TO INDEX %d\n",plumedgroup); 
  }
  else error->all("Illegal fix plumed command");
  // commit my special mpi type rvec
  MPI_Type_contiguous( 3 ,MPI_DOUBLE,&MPI_RVEC);
  MPI_Type_commit(&MPI_RVEC); 
  free(mss);
  free(chg);
  free(int_p);
//  fprintf(screen,"EXITING PLUMED PROC %d\n",me);   

//  MPI_Finalize();
//  abort();

  return;
}

/* ---------------------------------------------------------------------- */

FixPlumed::~FixPlumed()
{

}

/* ---------------------------------------------------------------------- */

int FixPlumed::setmask()
{
  // set with a bitmask how and when apply the force from plumed 
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPlumed::init()
{

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixPlumed::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixPlumed::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPlumed::post_force(int vflag)
{
  plumed_interface();
}

/* ---------------------------------------------------------------------- */

void FixPlumed::plumed_interface()
{
//printf("ENTERING PLUMED INTERFACE\n");
  int me,nprocs;
  MPI_Comm_rank(world,&me) ;
  MPI_Comm_size(world,&nprocs ) ;
 // printf("MY PROC IS %d SIZE IS %d\n",me,nprocs  );

  // alias the atoms

  double **f = atom->f;
  int   *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;

  // retrieve positions in a single vector
  int *my_atoms;
  int *my_idx,*my_backtable;
  rvec *my_pos,*my_force;
  int i,j;
  my_atoms=( int *)calloc(nprocs,sizeof(int));
  my_backtable=( int *)malloc(atom->natoms*sizeof(int));
  my_pos=( rvec *)malloc(atom->natoms*sizeof(rvec));
  my_force=( rvec *)malloc(atom->natoms*sizeof(rvec));
  my_idx=( int *)calloc(atom->natoms,sizeof(int));
  for (i = 0; i < nlocal; i++){ // loop over local atoms 
    if (mask[i] & plumedgroup2bit) {   // 
//      printf("LOCAL ATOM %d :I BELONG TO PLUMED IN PROC %d TAG: %d\n",i,me,tag[i]);
      // place it in the right vector:
      // tag to mypos :  index is needed
      j=my_atoms[me]; 
      my_pos[j][0]=x[i][0];  
      my_pos[j][1]=x[i][1];  
      my_pos[j][2]=x[i][2];  
      my_backtable[j]=i;
      my_idx[j]=tag[i]-1;  
 //     printf("PROC %d IDX %d LOCAL_IDX %d\n",me, my_idx[j],my_backtable[j]);
      my_atoms[me]++;
    }
  }
  // gather all the infos    
  int sum_allatoms, meatoms;
  meatoms=my_atoms[me];
  MPI_Allreduce(&meatoms,&sum_allatoms,1,MPI_INT,MPI_SUM,world);
  MPI_Allgather(&meatoms,1,MPI_INT,my_atoms,1,MPI_INT,world);
  // test the allgather 

//   for (i=0;i<nprocs;i++){
//      printf("----PROC %d RIGHT %d %d\n",me,i,my_atoms[i]);
//   }

  // put all the atoms in the same bucket: allocate it
//  printf("total number of atoms %d\n",sum_allatoms); 
  // incr pos contains the position in the buffer  
  int  *incr_pos;
  incr_pos=(int *)calloc(nprocs,sizeof(int));
  j=0;
  // only the root has informations
  if(me==0)for(i=0;i<nprocs;i++){incr_pos[i]=j;j+=my_atoms[i];}
  // this is the final vector to be allocated for store

  rvec *allpos, *allforce;
  int *allidx; 
  allpos=(rvec *)calloc(sum_allatoms,sizeof(rvec));
  allforce=(rvec *)calloc(sum_allatoms,sizeof(rvec));
  allidx=(int *)calloc(sum_allatoms,sizeof(int));
  MPI_Gatherv(my_pos,my_atoms[me],MPI_RVEC,allpos,&my_atoms[me],&incr_pos[me],MPI_RVEC,0,world); 
  MPI_Gatherv(my_idx,my_atoms[me],MPI_INT,allidx,&my_atoms[me],&incr_pos[me],MPI_INT,0,world); 


//if(me==0){
//   printf("PROC %d RIGHT \n",me);
//   for (i=0;i<nprocs;i++){
//      printf("P %d : ",i);
//      for (j=0;j<my_atoms[i];j++){
//            printf(" %d ",allidx[incr_pos[i]+j]);
//      }
//      printf("\n");
//   }
//   printf("\n");
//  }else{
////   printf("PROC %d WRONG \n",me);
////   for (i=0;i<nprocs;i++){
////      printf("P %d :  ",i);
////      for (j=0;j<my_atoms[i];j++){
////            printf(" %d ",allidx[incr_pos[i]+j]);
////      }
////      printf("\n");
////   }
////   printf("\n");
//  }
 

  
  if(!me){
   // for(i=0;i<sum_allatoms;i++){
   //     printf("COORD IDX %d COORD %f %f %f\n",allidx[i],allpos[i][0],allpos[i][1],allpos[i][2]);
   // }
    // now do the plumed calculation 
    plumed_ptr->meta_force_calculation(allidx,allpos,allforce,sum_allatoms,domain);

  } 
  // add the forces on each node 
  // scatter on the vector
  MPI_Scatterv(allforce,&my_atoms[me],&incr_pos[me],MPI_RVEC,my_force,my_atoms[me],MPI_RVEC,0,world); 
  // now each vector should have its force: loop on them
  for (i = 0; i < my_atoms[me] ; i++){ // loop over local atoms 
      // place it in the right vector:
      j=my_backtable[i];
      f[j][0]+=my_force[i][0];
      f[j][1]+=my_force[i][1];
      f[j][2]+=my_force[i][2];
  }

  // free all the vectors
  free(allpos); 
  free(allforce); 
  free(allidx); 
  free(my_atoms);
  free(my_pos);
  free(my_force);
  free(my_idx);
  free(my_backtable);
  free(incr_pos);
//  MPI_Finalize();
//  abort();


 
//  printf("EXITING  PLUMED INTERFACE\n");
};

/* ---------------------------------------------------------------------- */

void FixPlumed::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPlumed::min_post_force(int vflag)
{
  post_force(vflag);
}
EOF
cat >./fix_plumed.h  << \EOF
/*
*******************************************************************************
*                                                                             *
*                                PLUMED                                       *
*   A Portable Plugin for Free Energy Calculations with Molecular Dynamics    *
*                              VERSION 1.3                                    *
*                                                                             *
*******************************************************************************
*
*  
*  Copyright (c) 2008-2011 The PLUMED team.
*  See http://www.plumed-code.org for more information. 
*
*  This file is part of PLUMED.
*
*  PLUMED is free software: you can redistribute it and/or modify
*  it under the terms of the GNU Lesser General Public License as 
*  published by the Free Software Foundation, either version 3 of 
*  the License, or (at your option) any later version.
*
*  PLUMED is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General
*  Public License along with PLUMED.  
*  If not, see <http://www.gnu.org/licenses/>.
*
*  For more info, see:  http://www.plumed-code.org
*  or subscribe to plumed-users@googlegroups.com
*
*/

#ifdef FIX_CLASS

FixStyle(plumed,FixPlumed)

#else

#ifndef LMP_FIX_PLUMED_H
#define LMP_FIX_PLUMED_H

#include "fix.h"
// the plumed header that defines the class//
#include "metadyn.h"

namespace LAMMPS_NS {

class FixPlumed : public Fix {
 public:
  FixPlumed(class LAMMPS *, int, char **);
  ~FixPlumed();
  // class plumed included: a pointer that initialize the class and create all the interfaces 
  class Plumed *plumed;
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);

 private:
  Plumed *plumed_ptr; 
  int plumedgroup,plumedgroup2bit; 
  rvec *global_positions;
  int group2bit;
  int nlevels_respa;
  MPI_Datatype MPI_RVEC;
  void plumed_interface();
};

};

#endif
#endif

EOF

cd ../

#(this substitutes old "patch" command in diff directory)
mv Makefile Makefile.preplumed
sed 's/PACKUSER =/PACKUSER = user-plumed /' Makefile.preplumed > Makefile

for file in ./MAKE/Makefile.* ; do
  mv "$file" "${file}.preplumed"
  sed s/"CCFLAGS \="/"CCFLAGS \=  -DLAMMPS_PLUMED"/ ${file}.preplumed > $file
done
make yes-user-plumed
}

function to_do_before_revert () {
# create links
  echo > /dev/null
}

function to_do_after_revert () {
  make no-user-plumed
  rm -rf USER-PLUMED
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

