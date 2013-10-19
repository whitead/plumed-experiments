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

void PREFIX helix_restraint(int i_c, struct mtd_data_s *mtd_data) 
{
  real **psi;
  real tnorm[6], prdt[6], prdt_sign[6];
  real dpsi[6], dps[6], cos_psi[6];
  real nv1[6], nv2[6];
  real vp1[6][3], vp2[6][3];
  real r01[6][3], r21[6][3], r23[6][3];
  real psisum, psi_sum_tmp, pi;
  real ddotd[9], dnormd[9], dcosd[9];
  real fd1, fd2, fmod, mod_rij;
  int i, ii, jj, idtc, ndtc;
  int a0, a1, a2, a3;
 
  ndtc   = colvar.type[i_c];
  psisum = 0.;
  pi     = M_PI; 
  psi    = (real **)float_2d_array_alloc(6,ndtc);

  for(idtc=0;idtc<ndtc;idtc++) {

   psi_sum_tmp = 1.;

   for(ii=0;ii<6;ii++) {

    a0 = colvar.cvatoms[i_c][24*idtc+ii*4+0];
    a1 = colvar.cvatoms[i_c][24*idtc+ii*4+1];
    a2 = colvar.cvatoms[i_c][24*idtc+ii*4+2];
    a3 = colvar.cvatoms[i_c][24*idtc+ii*4+3];
   
    minimal_image(mtd_data->pos[a1],mtd_data->pos[a0],&mod_rij,r01[ii]);
    minimal_image(mtd_data->pos[a2],mtd_data->pos[a1],&mod_rij,r21[ii]);
    minimal_image(mtd_data->pos[a3],mtd_data->pos[a2],&mod_rij,r23[ii]);

    oprod(r01[ii], r21[ii], vp1[ii]);
    oprod(r21[ii], r23[ii], vp2[ii]);
    nv1[ii] = norm(vp1[ii]);
    nv2[ii] = norm(vp2[ii]);
    tnorm[ii] = nv1[ii] * nv2[ii];

    prdt[ii] = iprod(vp1[ii], vp2[ii]); 
    prdt_sign[ii] = iprod(r01[ii],vp2[ii]);
    prdt_sign[ii] /= fabs(prdt_sign[ii]);

    cos_psi[ii] = prdt[ii] / tnorm[ii];   
    if(fabs(cos_psi[ii]) > 1.) cos_psi[ii] /= fabs(cos_psi[ii]);

    psi[ii][idtc] = prdt_sign[ii] * acos(cos_psi[ii]);

    if(ii<3) dps[ii] = psi[ii][idtc]-colvar.map0[i_c][idtc];
    else     dps[ii] = psi[ii][idtc]-colvar.map1[i_c][idtc];

    psi_sum_tmp = psi_sum_tmp*0.5*(cos(dps[ii])+1.);

   } 

   psisum += psi_sum_tmp; 

   for(ii=0;ii<6;ii++) {     
    dpsi[ii]    = -0.5*sin(dps[ii]);
    for(jj=0;jj<6;jj++) if(ii!=jj) dpsi[ii] = dpsi[ii]*0.5*(cos(dps[jj])+1.); 
   }

   for(ii=0;ii<6;ii++) { 
    /* derivatives of the scalar product */
    ddotd[0] =  r21[ii][1]*vp2[ii][2] - r21[ii][2]*vp2[ii][1];
    ddotd[1] = -r21[ii][0]*vp2[ii][2] + r21[ii][2]*vp2[ii][0];
    ddotd[2] =  r21[ii][0]*vp2[ii][1] - r21[ii][1]*vp2[ii][0];

    ddotd[3] = -r21[ii][1]*vp1[ii][2] + r21[ii][2]*vp1[ii][1];
    ddotd[4] =  r21[ii][0]*vp1[ii][2] - r21[ii][2]*vp1[ii][0];
    ddotd[5] = -r21[ii][0]*vp1[ii][1] + r21[ii][1]*vp1[ii][0];

    ddotd[6] =  r23[ii][1]*vp1[ii][2] - r01[ii][1]*vp2[ii][2] - r23[ii][2]*vp1[ii][1] + r01[ii][2]*vp2[ii][1];
    ddotd[7] = -r23[ii][0]*vp1[ii][2] + r01[ii][0]*vp2[ii][2] + r23[ii][2]*vp1[ii][0] - r01[ii][2]*vp2[ii][0];
    ddotd[8] =  r23[ii][0]*vp1[ii][1] - r01[ii][0]*vp2[ii][1] - r23[ii][1]*vp1[ii][0] + r01[ii][1]*vp2[ii][0];

    /* derivatives of the norm part */
    fd1 = nv2[ii] / nv1[ii];
    fd2 = nv1[ii] / nv2[ii];

    dnormd[0] = -ddotd[3] * fd1;
    dnormd[1] = -ddotd[4] * fd1;
    dnormd[2] = -ddotd[5] * fd1;

    dnormd[3] = -ddotd[0] * fd2;
    dnormd[4] = -ddotd[1] * fd2;
    dnormd[5] = -ddotd[2] * fd2;

    dnormd[6] = ( r23[ii][1]*vp2[ii][2]-r23[ii][2]*vp2[ii][1] )*fd2 + ( r01[ii][2]*vp1[ii][1]-r01[ii][1]*vp1[ii][2] )*fd1;
    dnormd[7] = (-r23[ii][0]*vp2[ii][2]+r23[ii][2]*vp2[ii][0] )*fd2 + ( r01[ii][0]*vp1[ii][2]-r01[ii][2]*vp1[ii][0] )*fd1;
    dnormd[8] = ( r23[ii][0]*vp2[ii][1]-r23[ii][1]*vp2[ii][0] )*fd2 + (-r01[ii][0]*vp1[ii][1]+r01[ii][1]*vp1[ii][0] )*fd1;

    for(i=0;i<9;i++)
        dnormd[i] /= tnorm[ii]*tnorm[ii];

         /* derivatives of cosine */
    for(i=0;i<9;i++)
        dcosd[i] = ddotd[i]/tnorm[ii] - prdt[ii]*dnormd[i];

    /* check for numerical inconsistency */
    if(fabs(psi[ii][idtc])<0.00001 || fabs(psi[ii][idtc]-pi)<0.00001 || fabs(psi[ii][idtc]+pi)<0.00001)
        fmod = 0.;
    else
        fmod = prdt_sign[ii] / fabs(sin(psi[ii][idtc]));

    /* atom derivatives */
    colvar.myder[i_c][24*idtc+ii*4+0][0] = dcosd[0] * fmod * dpsi[ii];
    colvar.myder[i_c][24*idtc+ii*4+0][1] = dcosd[1] * fmod * dpsi[ii];
    colvar.myder[i_c][24*idtc+ii*4+0][2] = dcosd[2] * fmod * dpsi[ii];

    colvar.myder[i_c][24*idtc+ii*4+1][0] = -(dcosd[0]-dcosd[6]) * fmod * dpsi[ii];
    colvar.myder[i_c][24*idtc+ii*4+1][1] = -(dcosd[1]-dcosd[7]) * fmod * dpsi[ii];
    colvar.myder[i_c][24*idtc+ii*4+1][2] = -(dcosd[2]-dcosd[8]) * fmod * dpsi[ii];

    colvar.myder[i_c][24*idtc+ii*4+2][0] = -(dcosd[6]-dcosd[3]) * fmod * dpsi[ii];
    colvar.myder[i_c][24*idtc+ii*4+2][1] = -(dcosd[7]-dcosd[4]) * fmod * dpsi[ii];
    colvar.myder[i_c][24*idtc+ii*4+2][2] = -(dcosd[8]-dcosd[5]) * fmod * dpsi[ii];

    colvar.myder[i_c][24*idtc+ii*4+3][0] = -dcosd[3] * fmod * dpsi[ii];
    colvar.myder[i_c][24*idtc+ii*4+3][1] = -dcosd[4] * fmod * dpsi[ii];
    colvar.myder[i_c][24*idtc+ii*4+3][2] = -dcosd[5] * fmod * dpsi[ii];

   }
  }

  colvar.ss0[i_c] = psisum;

  free_2dr_array_alloc(psi,6);

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 int PREFIX read_helix   (char **word,int count,t_plumed_input *input,int *iline,FILE *fplog)
{
  int i, iat, iw;
  int atnr[24];
  int help;
  double delta, r_0, r_1;
  char string[400];

  help = 0;
  
  iw=seek_word(word,"NLOOP");
  if(iw>=0) {sscanf(word[iw+1],"%i", &colvar.type[count]);} else {fprintf(fplog,"|- NEEDED NLOOP KEYWORD FOR HELIX\n"); help=1;}
  iw=seek_word(word,"SIGMA");
  if(iw>=0) {
       sscanf(word[iw+1],"%lf", &delta);
       colvar.delta_r[count]  = (real) delta; }

   if(help){
         fprintf(fplog,"|- HELIX SYNTAX:\n");
         fprintf(fplog,"|- NLOOP              : number of loops\n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"e.g. \n");
         fprintf(fplog,"HELIX NLOOP 2 SIGMA 0.1 \n");
         fprintf(fplog,"12 14 20 22   22 24 35 37   37 39 49 51  -0.785    10 12 14 20  20 22 24 35   35 37 39 49   -1.200 \n");
         fprintf(fplog,"22 24 35 37   37 39 49 51   51 53 59 61  -0.785    20 22 24 35  35 37 39 49   49 51 53 59   -1.200 \n");
         plumed_error("PluMeD dead with errors: check log file");
  }


  colvar.natoms[count]   = 24*colvar.type[count]; 

  fprintf(fplog, "\n%1i-HELIX: %i LOOPS; ", count+1, colvar.type[count]);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n");

  snew(colvar.cvatoms[count], colvar.natoms[count]);
  snew(colvar.myder[count],   colvar.natoms[count]);
  snew(colvar.map0[count],    colvar.type[count]);
  snew(colvar.map1[count],    colvar.type[count]);

  for(iat=0;iat<colvar.type[count];iat++){
    (*iline)++;
    for(i=0;i<12;i++)  sscanf(input->words[*iline][i],"%i", &atnr[i]); 
    sscanf(input->words[*iline][12],"%lf",&r_0);
    for(i=13;i<25;i++) sscanf(input->words[*iline][i],"%i", &atnr[i-1]);
    sscanf(input->words[*iline][25],"%lf",&r_1);

    for(i=0;i<24;i++) colvar.cvatoms[count][24*iat+i]   = atnr[i]-1;
    colvar.map0[count][iat] = (real) r_0;
    colvar.map1[count][iat] = (real) r_1;

    fprintf(fplog, "|--LOOP %i, ATOMS: ",iat+1);
    for (i=0;i<3;i++) fprintf(fplog, "%i %i %i %i - ",atnr[i*4+0],atnr[i*4+1],atnr[i*4+2],atnr[i*4+3]);
    fprintf(fplog, " R_0 :: %lf\n",colvar.map0[count][iat]);
    fprintf(fplog, "          ATOMS: ");
    for (i=0;i<3;i++) fprintf(fplog, "%i %i %i %i - ",atnr[12+i*4+0],atnr[12+i*4+1],atnr[12+i*4+2],atnr[12+i*4+3]);
    fprintf(fplog, " R_1 :: %lf\n",colvar.map1[count][iat]); 
  }

  return colvar.type[count];

}
