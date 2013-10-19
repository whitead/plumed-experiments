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

void PREFIX dipole_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int  i, ai;
  rvec dipje;
  real dfunc;

  clear_rvec(dipje);

  for(i=0;i<colvar.natoms[i_c];i++) {
    ai = colvar.cvatoms[i_c][i];
    dipje[0] += mtd_data->pos[ai][0]*mtd_data->charge[ai];
    dipje[1] += mtd_data->pos[ai][1]*mtd_data->charge[ai];
    dipje[2] += mtd_data->pos[ai][2]*mtd_data->charge[ai];
  }

  colvar.ss0[i_c] = norm(dipje);

  for(i=0;i<colvar.natoms[i_c];i++) {
    ai = colvar.cvatoms[i_c][i];
    dfunc = mtd_data->charge[ai]/colvar.ss0[i_c];
    colvar.myder[i_c][i][0] = dfunc*dipje[0];
    colvar.myder[i_c][i][1] = dfunc*dipje[1];
    colvar.myder[i_c][i][2] = dfunc*dipje[2];
  }

}

//-------------------------------------------------------------------------------------------------------

int PREFIX read_dipole(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, j, iat, iw;
  double delta = 0.0;
  char string[400];
  int help;
  help = 0;

  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else {
    fprintf(fplog,"|- NEEDED LIST KEYWORD FOR DIPOLE\n");
    help=1;
  }

  iw=seek_word(word,"SIGMA");
 if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  if(help){
         fprintf(fplog,"|- DIPOLE SYNTAX:\n");
         fprintf(fplog,"|- LIST             : atom number/name of the group of atoms\n");
         fprintf(fplog,"|- SIGMA            : hills width for this cv\n");
         fprintf(fplog,"e.g. \n");
         fprintf(fplog,"DIPOLE LIST <g1> SIGMA 0.1 \n");
         fprintf(fplog, "         g1->    \n");
         fprintf(fplog, "         6 10 15 30 32    \n");
         fprintf(fplog, "         g1<-    \n"); 
         plumed_error("PluMeD dead with errors: check log file");
  } 
 
  colvar.type_s[count]   = 8;

  fprintf(fplog, "\n%1i-DIPOLE: ATOMS %i ", count+1, colvar.natoms[count]);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 

  fprintf(fplog,"|- SET MEMBERS: ");
  for(i=0;i<colvar.natoms[count];i++){fprintf(fplog," %d ",colvar.cvatoms[count][i]+1);if((i+1)%20==0)fprintf(fplog,"\n               ");}fprintf(fplog,"\n\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  return colvar.natoms[count];
}
