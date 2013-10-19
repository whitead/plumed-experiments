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

void PREFIX elstpot_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int firstAtom, secondAtom,iat, i, j, k, ix;
  real r_0=colvar.realpar[i_c][0][0];   
  real cut=colvar.realpar[i_c][1][0];   
  real elstpot;
  real mass1=0.;
  elstpot=0.; 

  for(i=0;i<colvar.natoms[i_c];i++){
    for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] = 0.;
  }

  real mod_rij,func,dfunc;
  rvec rij, sum1;
  // position of the point where it's felt the potential 
  sum1[0] = sum1[1] = sum1[2] = 0.;
  for(i=0;i<colvar.list[i_c][0];i++) {
    iat = colvar.cvatoms[i_c][i];
    sum1[0] += mtd_data->mass[iat]*mtd_data->pos[iat][0];
    sum1[1] += mtd_data->mass[iat]*mtd_data->pos[iat][1];
    sum1[2] += mtd_data->mass[iat]*mtd_data->pos[iat][2];
    mass1 += mtd_data->mass[iat];
  }
  sum1[0] /= mass1; sum1[1] /= mass1; sum1[2] /= mass1;
 
  for(j=colvar.list[i_c][0];j<colvar.natoms[i_c];j++) {
      iat = colvar.cvatoms[i_c][j];
      minimal_image(sum1, mtd_data->pos[iat], &mod_rij, rij);
      func=(mtd_data->charge[iat]/mod_rij);
      elstpot+=func*coscut(mod_rij,r_0,cut);
      //elstpot+=func;
      dfunc=-func*coscut(mod_rij,r_0,cut)/mod_rij+func*dcoscut(mod_rij,r_0,cut) ;
      //dfunc=-func/mod_rij ;
      for(ix=0;ix<3;ix++) {
           colvar.myder[i_c][j][ix] -= dfunc*rij[ix]/mod_rij;
      }
      for(i=0;i<colvar.list[i_c][0];i++) {
           k = colvar.cvatoms[i_c][i];
           for(ix=0;ix<3;ix++) {
        	 colvar.myder[i_c][i][ix] -= -dfunc*rij[ix]*mtd_data->mass[k]/(mod_rij*mass1);
           }
      }
  }
  colvar.ss0[i_c] = elstpot;

}

//--------------------------------------------------------------------------------------------------------------

int PREFIX read_elstpot(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw, iat, j, help;
  double r_0, cut,beta ;
  double delta = 0.0;
  char string[400];
  real threshold, value;

  help=0;

  iw = seek_word(word,"LIST");
  if(iw>=0){   
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;


  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR ELSTPOT\n"); help=1;}
  // the mu for fermi function 
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }
  // the beta for fermi function 
  iw=seek_word(word,"R_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &r_0); } else { fprintf(fplog,"|- NEEDED R_0 KEYWORD FOR ELSTPOT\n"); help=1;}
  // a clean cutoff
  iw=seek_word(word,"CUT");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &cut); }
  if(r_0>cut){fprintf(fplog,"|- NEEDED CUT > R_0 KEYWORD FOR ELSTPOT \n");}
  if(help){
          fprintf(fplog, "\n-ELSTPOT CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "ELSTPOT LIST <g1> <g2>  R_0 4.0  CUT 6.0  SIGMA 1.0 \n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "         5 1 6    \n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "       \n");
          fprintf(fplog, "         g2->    \n");
          fprintf(fplog, "         LOOP 6921 21786 5 \n");
          fprintf(fplog, "         g2<-    \n");
          fprintf(fplog, "       \n");
          plumed_error("PluMeD dead with errors: check log file");
  }


  (colvar.realpar[count][0][0])     = (real) r_0;
  (colvar.realpar[count][1][0])     = (real) cut;

  fprintf(fplog, "%1i-ELSTPOT ON  (1st SET: %i ATOMS) CAUSED BY (2nd SET: %i ATOMS); ", count+1, colvar.list[count][0], colvar.list[count][1]);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  fprintf(fplog, "|--FUNCTIONAL FORM:  sum_i  (q_2^i)/(r_2^i- sum_k( r_1^k ) ) f(r_2^i- sum_k( r_1^k ),r_0,cut ) Heaviside( r_2^i- sum_k( r_1^k )-cut)  \n");
  fprintf(fplog, "|--PARAMETERS: r_0= %f cut= %f \n",  colvar.realpar[count][0][0] ,colvar.realpar[count][1][0] );

  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  return colvar.natoms[count]; 
}
real PREFIX coscut(real r,real r_0, real cut){
   if(r<r_0)return 1.; 
   if(r>cut)return 0.;
   real d=(r-r_0)/(cut-r_0);
   return cos(0.5*M_PI*d); 
}
real PREFIX dcoscut(real r,real r_0, real cut){
   if(r<r_0)return 0.; 
   if(r>cut)return 0.;
   real d=(r-r_0)/(cut-r_0);
   return -sin(0.5*M_PI*d)*0.5*M_PI/(cut-r_0); 
}
