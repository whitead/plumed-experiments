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

void PREFIX torsion_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i, j, iat, k;
  real myp[4][3], totmasse[4], d[4][3];
  int t1, t2, t3, ix;
  rvec vp1, vp2, vp3;
  real pi, psi, nv1, nv2, n21, sc1, sc2;
  real twopi, in21, sign, cos_psi, sin_psi;
  rvec r01, r21, r23, m, n;
 
  twopi = M_2PI;
  pi = M_PI;

  totmasse[0] = totmasse[1] = totmasse[2] = totmasse[3] = 0.;
  myp[0][0] = myp[0][1] = myp[0][2] = 0.;
  myp[1][0] = myp[1][1] = myp[1][2] = 0.;
  myp[2][0] = myp[2][1] = myp[2][2] = 0.;
  myp[3][0] = myp[3][1] = myp[3][2] = 0.;
  k = 0;

  for(j=0;j<4;j++) {
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

#if defined (PLUMED_GROMACS)
#if defined (PLUMED_GROMACS45)
  psi = dih_angle(myp[0],myp[1],myp[2],myp[3],&mtd_data->metapbc,r01,r21,r23,m,n,
              &sign,&t1,&t2,&t3);
#else
  psi = dih_angle(myp[0],myp[1],myp[2],myp[3],&mtd_data->metapbc,r01,r21,r23,m,n,
              &cos_psi,&sign,&t1,&t2,&t3);
#endif
#else
  psi = dih_angle(myp[0],myp[1],myp[2],myp[3],r01,r21,r23,m,n,
              &cos_psi,&sign);
#endif

  oprod(r01,r21,vp1);
  oprod(r21,r23,vp2);
  oprod(vp1,vp2,vp3);

  if (psi >= pi) psi -= twopi;
  else if(psi < -pi) psi += twopi;

  nv1=norm(vp1);
  nv2=norm(vp2);
  n21=norm(r21);

  sc1=iprod(r01,r21);
  sc2=iprod(r23,r21);

  in21=1./(n21*n21);

  sc1 = sc1*in21;
  sc2 = sc2*in21;

  for(i=0;i<3;++i) {
    d[0][i]= -n21 * vp1[i] / (nv1*nv1);
    d[3][i]=  n21 * vp2[i] / (nv2*nv2);
    d[1][i]=  (sc1-1.0)*d[0][i] - sc2*d[3][i];
    d[2][i]=  (sc2-1.0)*d[3][i] - sc1*d[0][i];
  }

  // Now we do appropriate multiplication of the derivatives to get what we are interested in
  if( colvar.doTrig[i_c]==0 ){ colvar.ss0[i_c] = psi; }
  else if( colvar.doTrig[i_c]==1 ){
     cos_psi=cos( psi );
     for(i=0;i<3;++i){d[0][i]*= cos_psi; d[3][i]*= cos_psi; d[1][i]*= cos_psi; d[2][i]*= cos_psi;}
     colvar.ss0[i_c] = sin( psi );
  }
  else if( colvar.doTrig[i_c]==2 ){
     sin_psi=-sin( psi );
     for(i=0;i<3;++i){d[0][i]*= sin_psi; d[3][i]*= sin_psi; d[1][i]*= sin_psi; d[2][i]*= sin_psi;}
     colvar.ss0[i_c] = cos( psi );
  }
  else{ plumed_error("No trigonometric mode defined in torsion restraint"); }

  k=0;
  for(j=0;j<4;j++) {
    for(i=0;i<colvar.list[i_c][j];i++){
      iat = colvar.cvatoms[i_c][k];
      for(ix=0;ix<3;ix++) colvar.myder[i_c][k][ix] = -d[j][ix]*mtd_data->mass[iat]/totmasse[j];
      k++;
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------

int PREFIX read_torsion(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iat, j, iw, help;
  double delta = 0.0;
  char string[400];
  help=0;

  colvar.doTrig[count]=0;    // Default is to do no trigonometry (raw torsion)

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
             j=plumed_get_group(word[iw+4],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][3]=j;
  } else { fprintf(fplog,"|- NEEDED LIST KEYWORD FOR TORSION\n"); help=1;}
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

    if(help){
          fprintf(fplog, "\n-TORSION CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "TORSION LIST <g1> <g2> <g3> <g4> SIGMA 0.1\n");
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
          fprintf(fplog, "                 \n");
          fprintf(fplog, "         g4->    \n");
          fprintf(fplog, "         1 2     \n");
          fprintf(fplog, "         g3<-    \n"); 
          plumed_error("PluMeD dead with errors: check log file");          
  }
  iw=seek_word(word,"SIN");
  if(iw>=0){ colvar.doTrig[count]=1; }
  iw=seek_word(word,"COS");
  if(iw>=0){ colvar.doTrig[count]=2; }

  colvar.type_s[count]   = 5;

  if(colvar.doTrig[count]==0){
    fprintf(fplog, "\n%1i-TORSION: (1st SET: %i ATOMS), (2nd SET: %i ATOMS), (3rd SET: %i ATOMS) , (4th SET: %i ATOMS); ", count+1, 
      colvar.list[count][0], colvar.list[count][1], colvar.list[count][2], colvar.list[count][3]);
  }
  if(colvar.doTrig[count]==1){
    fprintf(fplog, "\n%1i-SINE OF TORSION: (1st SET: %i ATOMS), (2nd SET: %i ATOMS), (3rd SET: %i ATOMS) , (4th SET: %i ATOMS); ", count+1,
      colvar.list[count][0], colvar.list[count][1], colvar.list[count][2], colvar.list[count][3]);
  }
  if(colvar.doTrig[count]==2){
    fprintf(fplog, "\n%1i-COSINE OF TORSION: (1st SET: %i ATOMS), (2nd SET: %i ATOMS), (3rd SET: %i ATOMS) , (4th SET: %i ATOMS); ", count+1,
      colvar.list[count][0], colvar.list[count][1], colvar.list[count][2], colvar.list[count][3]);
  }
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 3rd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][2];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 4th SET MEMBERS: ");
  for(i=0;i<colvar.list[count][3];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  return colvar.natoms[count];
}
