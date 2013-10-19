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
#include "metadyn.h"
// Grant Rotskoff, 2013-06-06: grant.rotskoff@gmail.com //
/*
   LOCAL DENSITY COLLECTIVE VARIABLE

   Define a region in space with 4 atoms. The first 3 in the list form a plane.
   The 4th forms the diagonal of a cubiodal region. 

*/
// given a vector which is the position relative to the center
// of the box, calculate the number value of this position
// num_value = prod_{x,y,z} 
void PREFIX density_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  // the atoms defining the region
  int x0,x1,x2,x3;
  // the vectors defining the region
  // n is the normal to d0,d1
  rvec d0,d1,d2,n;
  // their moduli
  real mod_d0,mod_d1,mod_d2;
  // the cutoffs for the region
  real xc,yc,zc;
  // the box center
  rvec bcenter;

  int i;
  // initialize the derivatives
  for (i=0;i<colvar.natoms[i_c];i++) {
    colvar.myder[i_c][i][0] = 0.;
    colvar.myder[i_c][i][1] = 0.;
    colvar.myder[i_c][i][2] = 0.;
  }

  // get the atoms from the cv data structures
  x0=colvar.cvatoms[i_c][0];
  x1=colvar.cvatoms[i_c][1];
  x2=colvar.cvatoms[i_c][2];
  x3=colvar.cvatoms[i_c][3];
  
  
  // create the box region
  minimal_image(mtd_data->pos[x0],mtd_data->pos[x1],&mod_d0,d0);
  minimal_image(mtd_data->pos[x0],mtd_data->pos[x2],&mod_d1,d1);
  minimal_image(mtd_data->pos[x0],mtd_data->pos[x3],&mod_d2,d2);
  oprod(d0,d1,n);
  rvec xproj,yproj,zproj;
  xproj[0]=0.5*(iprod(d0,d2)/(mod_d0*mod_d0))*d0[0]; // project d2 onto d0
  xproj[1]=0.5*(iprod(d0,d2)/(mod_d0*mod_d0))*d0[1]; 
  xproj[2]=0.5*(iprod(d0,d2)/(mod_d0*mod_d0))*d0[2]; 
  yproj[0]=0.5*(iprod(d1,d2)/(mod_d1*mod_d1))*d1[0]; // project d2 onto d1
  yproj[1]=0.5*(iprod(d1,d2)/(mod_d1*mod_d1))*d1[1]; 
  yproj[2]=0.5*(iprod(d1,d2)/(mod_d1*mod_d1))*d1[2]; 
  zproj[0]=0.5*(iprod(n,d2)/(norm(n))*n[0]); // project d2 onto n
  zproj[1]=0.5*(iprod(n,d2)/(norm(n))*n[1]); 
  zproj[2]=0.5*(iprod(n,d2)/(norm(n))*n[2]); 
  bcenter[0]=mtd_data->pos[x0][0]+norm(xproj);
  bcenter[1]=mtd_data->pos[x0][1]+norm(yproj);
  bcenter[2]=mtd_data->pos[x0][2]+norm(zproj);
  //bcenter[1]=x0[1]+norm(yproj);
  //bcenter[2]=x0[2]+norm(zproj);
  // same for x_cutoff and y_cutoff
  xc=fabs(mtd_data->pos[x0][0]-bcenter[0]); // distance to from the center to x0 in the xdirection
  yc=fabs(mtd_data->pos[x0][1]-bcenter[1]); // distance to from the center to x0 in the xdirection
  zc=fabs(mtd_data->pos[x0][2]-bcenter[2]); // distance to from the center to x0 in the xdirection
  // get the z_cutoff by getting the distance to the plane defined by d0,d1
  // we project 1/2*d onto the normal vector to the plane
  // first, compute the cross product of d0, d1

  // the waters are the atoms after the first 4 atoms in the CV list
  int k,awater;
  real num_waters=0.;
  rvec d; // position of the water relative to the box center
  real fx,fy,fz,dfx,dfy,dfz,mod_d; // values of the number density
  for (k=colvar.list[i_c][0];k<colvar.natoms[i_c];k++) {
    awater=colvar.cvatoms[i_c][k];
    // check if the water lives within the region
    minimal_image(mtd_data->pos[awater], bcenter, &mod_d, d);
    // if there's a water, compute it's contribution to the density
    if (fabs(d[0])<xc && fabs(d[1])<yc && fabs(d[2])<zc) {
      fx=(1-(d[0]/xc)*(d[0]/xc))*(1-(d[0]/xc)*(d[0]/xc));
      fy=(1-(d[1]/yc)*(d[1]/yc))*(1-(d[1]/yc)*(d[1]/yc));
      fz=(1-(d[2]/zc)*(d[2]/zc))*(1-(d[2]/zc)*(d[2]/zc));
      num_waters+=(fx*fy*fz);
      // compute the derivatives to get the forces
      dfx=-4*fy*fz*(d[0]/xc-(d[0]/xc)*(d[0]/xc));
      dfy=-4*fx*fz*(d[0]/yc-(d[0]/yc)*(d[0]/yc));
      dfz=-4*fx*fy*(d[0]/zc-(d[0]/zc)*(d[0]/zc));
      colvar.myder[i_c][k][0]+=dfx;
      colvar.myder[i_c][k][1]+=dfy;
      colvar.myder[i_c][k][2]+=dfz;
    }
  }
  // update the value of the density
  colvar.ss0[i_c] = num_waters;
}
    
int PREFIX read_density(char **word, int cvnum, t_plumed_input *input, FILE *log)
{
  int iw,i,j,k;
  double delta=0.;
  int help=0;
  
  // get the region defining atom list
  iw=seek_word(word,"LIST");
  if (iw>=0) {
    // list 1 is the four atoms defining the region
    j=plumed_get_group(word[iw+1],&colvar.cvatoms[cvnum],colvar.natoms[cvnum],input,log);
    colvar.natoms[cvnum]+=j;
    colvar.list[cvnum][0]+=j;
    // list 2 is the solvent
    j=plumed_get_group(word[iw+2],&colvar.cvatoms[cvnum],colvar.natoms[cvnum],input,log);
    colvar.natoms[cvnum]+=j;
    colvar.list[cvnum][1]+=j;
  } else { fprintf(log,"|- DENSITY COLVAR REQUIRES A LIST OF 4 ATOMS FOLLOWED BY SOLVENT\n"); help=1;}
  
  // get the hill width
  iw=seek_word(word,"SIGMA");
  if (iw>=0) { 
    sscanf(word[iw+1],"%lf",&delta);
    colvar.delta_r[cvnum]=(real) delta;
  }

  if(help){
    fprintf(log,"|- DENSITY SYNTAX:\n");
    fprintf(log,"|- LIST               : atoms defining the region\n");
    fprintf(log,"|-                    : first 3 define a plane\n");
    fprintf(log,"|-                    : first and last define box diagonal\n");
    fprintf(log,"|-\n");
    fprintf(log,"|- LOOP               : water oxygens in the system\n");
    fprintf(log,"e.g. \n");
    fprintf(log,"DENSITY LIST <region> <solvent> SIGMA 0.1\n");
    fprintf(log, "         region->    \n");
    fprintf(log, "         17 29 88 100  \n");
    fprintf(log, "         region<-    \n");
    fprintf(log, "                 \n");
    fprintf(log, "         solvent->    \n");
    fprintf(log, "         LOOP 100 1000 3\n");
    fprintf(log, "         solvent<-    \n");
    plumed_error("PluMeD dead with errors: check log file");
  }

  colvar.type_s[cvnum] = 21;
  
  fprintf(log,"\n|- CV-%i DENSITY: (REGION DEFINED BY LIST 0: %d); ",cvnum+1,colvar.list[cvnum][0]);
  if (logical.do_hills) fprintf(log," SIGMA %f\n",colvar.delta_r[cvnum]);
  else fprintf(log,"\n");

  snew(colvar.myder[cvnum], colvar.natoms[cvnum]);
  // check that a selected water isn't defining the region
  for (k=0;k<colvar.list[cvnum][1];k++) {
    for (i=0;i<colvar.list[cvnum][0];i++) {
      if (colvar.cvatoms[cvnum][i]==colvar.cvatoms[cvnum][colvar.list[cvnum][0]+k]) {
        char buf[1024];
        sprintf(buf,"DENSITY: Atom %d defines the region. It may not be included in the solvent.",
                colvar.cvatoms[cvnum][i]);
        plumed_error(buf);
      }
    }
  }
  return colvar.list[cvnum][0];
}
