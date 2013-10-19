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

void PREFIX angle_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i, j, iat, k, ix;
  real myp[3][3], totmasse[3], d[3][3];
  rvec rij, sij, tij;
  real mod_rij, mod_sij;
  real eps, sign0;
  real caa, cbb, cab, ccc, sab;
  real ac, Vac, bc, dVac, cos_psi, sin_psi;

  totmasse[0] = totmasse[1] = totmasse[2] = 0.;
  myp[0][0] = myp[0][1] = myp[0][2] = 0.;
  myp[1][0] = myp[1][1] = myp[1][2] = 0.;
  myp[2][0] = myp[2][1] = myp[2][2] = 0.;
  k = 0;

  for(j=0;j<3;j++) {
    for(i=0;i<colvar.list[i_c][j];i++){
      iat = colvar.cvatoms[i_c][k];
      myp[j][0] += mtd_data->mass[iat]*mtd_data->pos[iat][0];
      myp[j][1] += mtd_data->mass[iat]*mtd_data->pos[iat][1];
      myp[j][2] += mtd_data->mass[iat]*mtd_data->pos[iat][2];
      totmasse[j] += mtd_data->mass[iat];
      k++;
    }
    myp[j][0] /= totmasse[j];
    myp[j][1] /= totmasse[j];
    myp[j][2] /= totmasse[j];
  }

  minimal_image(myp[1], myp[0], &mod_rij, rij);
  minimal_image(myp[1], myp[2], &mod_sij, sij);
  
  caa = norm2(rij);
  cab = iprod(rij,sij);
  cbb = norm2(sij);

  oprod(rij,sij,tij);
  sab = norm(tij);
  ccc = 1.0/sqrt(caa*cbb);
  ac = cab*ccc; // cosine of teta
  Vac = acos(ac);
  dVac = -ccc/sqrt(1.-ac*ac);

  for(ix=0;ix<3;ix++) {
    d[0][ix] = -dVac*(-cab/caa*rij[ix]+sij[ix]);
    d[1][ix] = dVac*(-cab/caa*rij[ix]-cab/cbb*sij[ix]+rij[ix]+sij[ix]);
    d[2][ix] = dVac*(cab/cbb*sij[ix]-rij[ix]);
  }

  // Now we do appropriate multiplication of the derivatives to get what we are interested in
  if( colvar.doTrig[i_c]==0 ){ colvar.ss0[i_c] = Vac; }
  else if( colvar.doTrig[i_c]==1 ){
     cos_psi=cos( Vac );
     for(i=0;i<3;++i){d[0][i]*= cos_psi;  d[1][i]*= cos_psi; d[2][i]*= cos_psi;}
     colvar.ss0[i_c] = sin( Vac );
  }
  else if( colvar.doTrig[i_c]==2 ){
     sin_psi=-sin( Vac );
     for(i=0;i<3;++i){d[0][i]*= sin_psi;  d[1][i]*= sin_psi; d[2][i]*= sin_psi;}
     colvar.ss0[i_c] = cos( Vac );
  }
  else{ plumed_error("No trigonometric mode defined in torsion restraint"); }



  k=0;
  for(j=0;j<3;j++) {
    for(i=0;i<colvar.list[i_c][j];i++){
      iat = colvar.cvatoms[i_c][k];
      for(ix=0;ix<3;ix++) colvar.myder[i_c][k][ix] = d[j][ix]*mtd_data->mass[iat]/totmasse[j];
      k++;
    }
  }

  //colvar.ss0[i_c] = Vac;
}

//---------------------------------------------------------------------------------------------------------------

int PREFIX read_angle(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i,j, iat, iw, help;
  double delta = 0.0;
  char string[400];
  
  help = 0;

  colvar.doTrig[count]=0;    // Default is to do no trigonometry (raw angle)

  iw = seek_word(word,"LIST");
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
  } else { fprintf(fplog,"|- NEEDED LIST KEYWORD FOR ANGLE\n"); help=1;}

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);  
             colvar.delta_r[count]  = (real) delta; }
 
  if(help){
          fprintf(fplog, "\n-ANGLE CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "      ANGLE LIST 12 17 20 SIGMA 1.0 \n");
          fprintf(fplog, "              \n");
          fprintf(fplog, "or in case of group: \n");
          fprintf(fplog, "      ANGLE LIST <g1> <g2> <g3> SIGMA 1.0  \n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "         6 10    \n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "                 \n");
          fprintf(fplog, "         g2->    \n");
          fprintf(fplog, "         8 15 21 \n");
          fprintf(fplog, "         g2<-    \n");
          fprintf(fplog, "                 \n"); 
          fprintf(fplog, "         g3->    \n");
          fprintf(fplog, "         23 29 31\n");
          fprintf(fplog, "         g3<-    \n");
          fprintf(fplog, "   \n");
          plumed_error("PluMeD dead with errors: check log file");
  }

  iw=seek_word(word,"SIN");
  if(iw>=0){ colvar.doTrig[count]=1; }
  iw=seek_word(word,"COS");
  if(iw>=0){ colvar.doTrig[count]=2; }

  colvar.type_s[count]   = 4;

  if(colvar.doTrig[count]==0){
	  fprintf(fplog, "\n%1i-ANGLE: (1st SET: %i ATOMS), (2nd SET: %i ATOMS), (3rd SET: %i ATOMS); ", count+1,
	    colvar.list[count][0], colvar.list[count][1], colvar.list[count][2]);
  }

  if(colvar.doTrig[count]==1){
	  fprintf(fplog, "\n%1i-SINE OF ANGLE: (1st SET: %i ATOMS), (2nd SET: %i ATOMS), (3rd SET: %i ATOMS); ", count+1,
	    colvar.list[count][0], colvar.list[count][1], colvar.list[count][2]);
  }

  if(colvar.doTrig[count]==1){
	  fprintf(fplog, "\n%1i-COSINE OF ANGLE: (1st SET: %i ATOMS), (2nd SET: %i ATOMS), (3rd SET: %i ATOMS); ", count+1,
	    colvar.list[count][0], colvar.list[count][1], colvar.list[count][2]);
  }

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

  return colvar.natoms[count];
}
