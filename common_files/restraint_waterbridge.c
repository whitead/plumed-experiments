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
//-----------------------------------------------------------------------------
// COLVAR = number of water molecules bridging between two groups of atoms
// \sum_w \sum_a \sum_b n_aw*n_bw   (w in OW, a in group A, b in group B)
//-----------------------------------------------------------------------------

#include "metadyn.h"

void PREFIX waterbridge_restraint(int i_c, struct mtd_data_s *mtd_data) 
{
  int firstAtom, secondAtom, middleAtom, i, j, k, ix;
  rvec rik, rjk;
  real mod_rik, mod_rjk, s_rik, s_rjk;
  real num, den, func1, dfunc1, func2, dfunc2;
  real r_0 = colvar.r_0[i_c];
  real r_safe;
  int nn = colvar.nn[i_c];
  int mm = colvar.mm[i_c];
  int bdoneik;
  real ncoord = 0., ncoord_tmp;

  func1 = 0; func2 = 0; dfunc1 = 0; dfunc2 = 0;
  r_safe = r_0*4.; 

  for(i=0;i<colvar.natoms[i_c];i++) {
    colvar.myder[i_c][i][0] = 0.;
    colvar.myder[i_c][i][1] = 0.;
    colvar.myder[i_c][i][2] = 0.;
  }  

  // k identifies water atoms
  for(k=colvar.list[i_c][0]+colvar.list[i_c][1];k<colvar.natoms[i_c];k++){
    middleAtom = colvar.cvatoms[i_c][k];
    ncoord_tmp = 0.;
    for (i=0;i<colvar.list[i_c][0];i++) {
      firstAtom = colvar.cvatoms[i_c][i];
      // n_ik and gradients are not yet calculated
      bdoneik = 0;
      minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[middleAtom], &mod_rik, rik);
      // calculation continues only if i is near to k
      if(mod_rik < r_safe) {
        for(j=colvar.list[i_c][0];j<colvar.list[i_c][0]+colvar.list[i_c][1];j++) {
          secondAtom = colvar.cvatoms[i_c][j];
          minimal_image(mtd_data->pos[secondAtom], mtd_data->pos[middleAtom], &mod_rjk, rjk);
          // calculation continues only if j is near to k
          if(mod_rjk < r_safe) {
            // n_ik and gradients are calculated if undone before
            if(bdoneik==0) {
              s_rik = mod_rik/r_0;
              num = (1.-pow(s_rik, nn));
              den = (1.-pow(s_rik, mm));
              func1 = num/den;
              dfunc1 = (-nn*pow(s_rik,nn-1)/den+num/den/den*mm*pow(s_rik,mm-1))/mod_rik/r_0;
              bdoneik = 1;
            }
            // n_jk and gradients are calculated
            s_rjk = mod_rjk/r_0;
            num = (1.-pow(s_rjk, nn));
            den = (1.-pow(s_rjk, mm));
            func2 = num/den;
            dfunc2 = (-nn*pow(s_rjk,nn-1)/den+num/den/den*mm*pow(s_rjk,mm-1))/mod_rjk/r_0;

            // final result
            ncoord += func1*func2;
            ncoord_tmp += func1*func2;

            // dfunc : force (-derivative)
            for(ix=0;ix<3;ix++) { 
              colvar.myder[i_c][i][ix] += dfunc1*rik[ix]*func2;
              colvar.myder[i_c][j][ix] += dfunc2*rjk[ix]*func1;
              colvar.myder[i_c][k][ix] -= dfunc1*rik[ix]*func2 + dfunc2*rjk[ix]*func1;
            }
          }
        }
      }
    }
  }

  colvar.ss0[i_c] = ncoord;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int PREFIX read_waterbridge(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iat, iw, j, k;
  double delta = 0.0;
  double r_0;
  char string[400];
  int help;
  help = 0;


  iw=seek_word(word,"LIST");
  if(iw>=0) { 
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;
             j=plumed_get_group(word[iw+3],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][2]=j;
  } else { fprintf(fplog,"|- NEEDED LIST KEYWORD FOR WATERBRIDGE\n"); help=1;}
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }
  iw=seek_word(word,"NN");
  if(iw>=0) { 
      sscanf(word[iw+1],"%i", &colvar.nn[count]);
  } else {
    fprintf(fplog,"|- NEEDED NN KEYWORD FOR WATERBRIDGE\n");
    help=1;
  }
  iw=seek_word(word,"MM");
  if(iw>=0) { 
      sscanf(word[iw+1],"%i", &colvar.mm[count]);
  } else {
    fprintf(fplog,"|- NEEDED MM KEYWORD FOR WATERBRIDGE\n");
    help=1;
  }
  iw=seek_word(word,"R_0");
  if(iw>=0) { 
     sscanf(word[iw+1],"%lf", &r_0);
  } else {
    fprintf(fplog,"|- NEEDED R_0 KEYWORD FOR WATERBRIDGE\n");
    help=1;
  }

   if(help){
         fprintf(fplog,"|- WATERBRIDGE SYNTAX:\n");
         fprintf(fplog,"|- LIST               : groups name\n");
         fprintf(fplog,"|- NN                 : switching function parameter\n");
         fprintf(fplog,"|- MM                 : switching function parameter\n");
         fprintf(fplog,"|- R_0                : switching function parameter\n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"e.g. \n");
         fprintf(fplog,"WATERBRIDGE LIST <type1> <type2> <solvent> NN 8 MM 12 R_0 4.0 SIGMA 0.1 \n");
         fprintf(fplog, "         type1->    \n");
         fprintf(fplog, "         6 10    \n");
         fprintf(fplog, "         type1<-    \n");
         fprintf(fplog, "                 \n");
         fprintf(fplog, "         type2->    \n");
         fprintf(fplog, "         8 15 21 \n");
         fprintf(fplog, "         type2<-    \n");
         fprintf(fplog, "                 \n");
         fprintf(fplog, "         solvent->    \n");
         fprintf(fplog, "         LOOP 100 1000 3\n");
         fprintf(fplog, "         solvent<-    \n");
         plumed_error("PluMeD dead with errors: check log file");
  }

  colvar.type_s[count] = 20;
  colvar.r_0[count]      = (real) r_0; 
  
  fprintf(fplog, "\n%1i-WATERBRIDGE: (1st SET: %i ATOMS), (2nd SET: %i ATOMS), (3rd SET: %i ATOMS); ", count+1,
    colvar.list[count][0], colvar.list[count][1], colvar.list[count][2]);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 3rd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][2];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  // checking for repeated atoms
  for(k=0;k<colvar.list[count][2];k++) {
    for(i=0;i<colvar.list[count][0]+colvar.list[count][1];i++) {
      if(colvar.cvatoms[count][i]==colvar.cvatoms[count][colvar.list[count][0]+colvar.list[count][1]+k]) {
       char buf[1024];
       sprintf(buf, "WATERBRIDGE: atom %d is a water oxygen",colvar.cvatoms[count][i]);
       plumed_error(buf); 
      }
    }
  }

  return colvar.list[count][0]+colvar.list[count][1];
}
