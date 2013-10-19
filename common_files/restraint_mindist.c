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

void PREFIX mindist_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int firstAtom, secondAtom, i, j, ix;
  rvec rij;
  real mod_rij;
  double sum, pe, beta, dist;
  sum=0.0;

  beta = colvar.beta[i_c];

  for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] = 0.;

  for(i=0;i<colvar.list[i_c][0];i++){
    firstAtom  = colvar.cvatoms[i_c][i];
    for(j=colvar.list[i_c][0];j<colvar.natoms[i_c];j++){
      secondAtom = colvar.cvatoms[i_c][j];
      if(firstAtom != secondAtom){
        if(colvar.cell_pbc[i_c]){
          minimal_image(mtd_data->pos[firstAtom],mtd_data->pos[secondAtom], &mod_rij, rij);
        } else {
          rij[0] = mtd_data->pos[firstAtom][0] - mtd_data->pos[secondAtom][0];
          rij[1] = mtd_data->pos[firstAtom][1] - mtd_data->pos[secondAtom][1];
          rij[2] = mtd_data->pos[firstAtom][2] - mtd_data->pos[secondAtom][2];
          mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
        };
        pe = exp(beta/mod_rij);
        sum += pe;
        for(ix=0;ix<3;ix++) { 
          colvar.myder[i_c][j][ix] += -pe*rij[ix]/(mod_rij*mod_rij*mod_rij);
          colvar.myder[i_c][i][ix] += pe*rij[ix]/(mod_rij*mod_rij*mod_rij);
        }
      }
    }
  }

  dist = beta/log(sum);
  colvar.ss0[i_c] = (real) dist;

  for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] = colvar.myder[i_c][i][ix]*dist*dist/sum;
}

//----------------------------------------------------------------------------------------------------------------

int PREFIX read_mindist(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw, j, iat, help;
  double delta = 0.0;
  char string[400];

  if(sizeof(real)==sizeof(float)) colvar.beta[count] = 20.;
  else colvar.beta[count] = 50.;

  help=0;

  colvar.cell_pbc[count]=1; // default is PBC

  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}
  iw=seek_word(word,"BETA");
  if(iw>=0){sscanf(word[iw+1],"%lf", &delta);
             colvar.beta[count]  = (real) delta; }

  iw = seek_word(word,"LIST");
  if(iw>=0){   
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR DISTANCE\n"); help=1;}
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  if(help){
          fprintf(fplog, "\n-MINDIST CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "MINDIST LIST <g1> <g2>  SIGMA 0.1 BETA 500.\n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "         6 10    \n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "                 \n");
          fprintf(fplog, "         g2->    \n");
          fprintf(fplog, "         8 15 21 \n");
          fprintf(fplog, "         g2<-    \n");
          fprintf(fplog, "                 \n"); 
          plumed_error("PluMeD dead with errors: check log file");
  }

  colvar.type_s[count] = 2;

  fprintf(fplog, "\n%1i-MINDIST: (1st SET: %i ATOMS), (2nd SET: %i ATOMS); ", count+1, colvar.list[count][0], colvar.list[count][1]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON; ");
  else                       fprintf(fplog, " PBC OFF; ");
  fprintf(fplog," BETA %f; ",colvar.beta[count]);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n");
  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  return colvar.natoms[count];
}
