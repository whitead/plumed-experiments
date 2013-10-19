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

void PREFIX alfabeta_restraint(int i_c, struct mtd_data_s *mtd_data) 
{
  real *psi;
  real psisum, dpsi, dps, cos_psi;
  real nv1, nv2, mod_rij;
  int i,ix, idtc, ndtc;
  int a0, a1, a2, a3;
  rvec vp1, vp2;
  rvec r01, r21, r23;
  real tnorm, prdt,prdt_sign,pi;
  real ddotd[9];
  real dnormd[9];
  real dcosd[9];
  real fd1, fd2, fmod;

  ndtc = colvar.type[i_c];
  snew(psi, ndtc);
  psisum = 0;
  pi = acos(-1.0);

  for(idtc=0;idtc<ndtc;idtc++) {

    a0 = colvar.cvatoms[i_c][4*idtc+0];
    a1 = colvar.cvatoms[i_c][4*idtc+1];
    a2 = colvar.cvatoms[i_c][4*idtc+2];
    a3 = colvar.cvatoms[i_c][4*idtc+3];

    minimal_image(mtd_data->pos[a1],mtd_data->pos[a0],&mod_rij,r01);
    minimal_image(mtd_data->pos[a2],mtd_data->pos[a1],&mod_rij,r21);
    minimal_image(mtd_data->pos[a3],mtd_data->pos[a2],&mod_rij,r23);

    oprod(r01, r21, vp1);
    oprod(r21, r23, vp2);
    nv1 = norm(vp1);
    nv2 = norm(vp2);
    tnorm = nv1 * nv2;

    prdt = iprod(vp1, vp2); 
    prdt_sign = iprod(r01,vp2);
    prdt_sign /= fabs(prdt_sign);

    cos_psi = prdt / tnorm;   
    if(fabs(cos_psi) > 1.) cos_psi /= fabs(cos_psi);

    psi[idtc] = prdt_sign * acos(cos_psi);
    dps = psi[idtc]-colvar.map0[i_c][idtc];
    psisum += 0.5*(cos(dps)+1.);
    dpsi = -0.5*sin(dps);

    /* derivatives of the scalar product */
    ddotd[0] =  r21[1]*vp2[2] - r21[2]*vp2[1];
    ddotd[1] = -r21[0]*vp2[2] + r21[2]*vp2[0];
    ddotd[2] =  r21[0]*vp2[1] - r21[1]*vp2[0];

    ddotd[3] = -r21[1]*vp1[2] + r21[2]*vp1[1];
    ddotd[4] =  r21[0]*vp1[2] - r21[2]*vp1[0];
    ddotd[5] = -r21[0]*vp1[1] + r21[1]*vp1[0];

    ddotd[6] =  r23[1]*vp1[2] - r01[1]*vp2[2] - r23[2]*vp1[1] + r01[2]*vp2[1];
    ddotd[7] = -r23[0]*vp1[2] + r01[0]*vp2[2] + r23[2]*vp1[0] - r01[2]*vp2[0];
    ddotd[8] =  r23[0]*vp1[1] - r01[0]*vp2[1] - r23[1]*vp1[0] + r01[1]*vp2[0];

    /* derivatives of the norm part */
    fd1 = nv2 / nv1;
    fd2 = nv1 / nv2;

    dnormd[0] = -ddotd[3] * fd1;
    dnormd[1] = -ddotd[4] * fd1;
    dnormd[2] = -ddotd[5] * fd1;

    dnormd[3] = -ddotd[0] * fd2;
    dnormd[4] = -ddotd[1] * fd2;
    dnormd[5] = -ddotd[2] * fd2;

    dnormd[6] = ( r23[1]*vp2[2]-r23[2]*vp2[1] )*fd2 + ( r01[2]*vp1[1]-r01[1]*vp1[2] )*fd1;
    dnormd[7] = (-r23[0]*vp2[2]+r23[2]*vp2[0] )*fd2 + ( r01[0]*vp1[2]-r01[2]*vp1[0] )*fd1;
    dnormd[8] = ( r23[0]*vp2[1]-r23[1]*vp2[0] )*fd2 + (-r01[0]*vp1[1]+r01[1]*vp1[0] )*fd1;

    for(i=0;i<9;i++)
        dnormd[i] /= tnorm*tnorm;

         /* derivatives of cosine */
    for(i=0;i<9;i++)
        dcosd[i] = ddotd[i]/tnorm - prdt*dnormd[i];

    /* check for numerical inconsistency */
    if(fabs(psi[idtc])<0.00001 || fabs(psi[idtc]-pi)<0.00001 || fabs(psi[idtc]+pi)<0.00001)
        fmod = 0.;
    else
        fmod = prdt_sign / fabs(sin(psi[idtc]));

    /* atom derivatives */
    colvar.myder[i_c][4*idtc+0][0] = dcosd[0] * fmod * dpsi;
    colvar.myder[i_c][4*idtc+0][1] = dcosd[1] * fmod * dpsi;
    colvar.myder[i_c][4*idtc+0][2] = dcosd[2] * fmod * dpsi;

    colvar.myder[i_c][4*idtc+1][0] = -(dcosd[0]-dcosd[6]) * fmod * dpsi;
    colvar.myder[i_c][4*idtc+1][1] = -(dcosd[1]-dcosd[7]) * fmod * dpsi;
    colvar.myder[i_c][4*idtc+1][2] = -(dcosd[2]-dcosd[8]) * fmod * dpsi;

    colvar.myder[i_c][4*idtc+2][0] = -(dcosd[6]-dcosd[3]) * fmod * dpsi;
    colvar.myder[i_c][4*idtc+2][1] = -(dcosd[7]-dcosd[4]) * fmod * dpsi;
    colvar.myder[i_c][4*idtc+2][2] = -(dcosd[8]-dcosd[5]) * fmod * dpsi;

    colvar.myder[i_c][4*idtc+3][0] = -dcosd[3] * fmod * dpsi;
    colvar.myder[i_c][4*idtc+3][1] = -dcosd[4] * fmod * dpsi;
    colvar.myder[i_c][4*idtc+3][2] = -dcosd[5] * fmod * dpsi;

  }

  colvar.ss0[i_c] = psisum;


  sfree(psi);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int PREFIX read_alfabeta(char **word, int count, t_plumed_input *input,int *iline,FILE *fplog)
{
  int iat, atnr1, atnr2, atnr3, atnr4, iw;
  int help;
  double delta = 0.0;
  double r_0;
  char string[400];

  help = 0;
  
  iw=seek_word(word,"NDIH");
  if(iw>=0) {
       sscanf(word[iw+1],"%i", &colvar.type[count]);
  } else {
    fprintf(fplog,"|- NEEDED NDIH KEYWORD FOR ALPHABETA\n");
    help=1;
  }
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

   if(help){
         fprintf(fplog,"|- ALPHABETA SYNTAX:\n");
         fprintf(fplog,"|- NDIH               : number of dihedrals\n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"e.g. \n");
         fprintf(fplog,"ALPHABETA NDIH 1 SIGMA 0.1 \n");
         fprintf(fplog,"168 170 172 188 3.14\n");
         plumed_error("PluMeD dead with errors: check log file"); 
  }

  colvar.type_s[count]   = 6;
  colvar.natoms[count]   = 4*colvar.type[count]; 

  fprintf(fplog, "\n%1i-ALPHABETA: %i DIHEDRAL; ", count+1, colvar.type[count]);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n");

  snew(colvar.cvatoms[count], colvar.natoms[count]);
  snew(colvar.myder[count],   colvar.natoms[count]);
  snew(colvar.map0[count],    colvar.type[count]);

  for(iat=0;iat<colvar.type[count];iat++){
    (*iline)++;
    sscanf(input->words[*iline][0],"%i",&atnr1);
    sscanf(input->words[*iline][1],"%i",&atnr2);
    sscanf(input->words[*iline][2],"%i",&atnr3);
    sscanf(input->words[*iline][3],"%i",&atnr4);
    sscanf(input->words[*iline][4],"%lf",&r_0);
    atnr1--;
    atnr2--;
    atnr3--;
    atnr4--;
    colvar.cvatoms[count][4*iat] = atnr1;
    colvar.cvatoms[count][4*iat+1] = atnr2;
    colvar.cvatoms[count][4*iat+2] = atnr3;
    colvar.cvatoms[count][4*iat+3] = atnr4;
    colvar.map0[count][iat] = (real) r_0;
    fprintf(fplog, "|--DIH %i, ATOMS: %i %i %i %i  R_0: %lf \n", iat+1, atnr1+1, atnr2+1, atnr3+1, atnr4+1, r_0);
  }
  fprintf(fplog,"\n");
  return colvar.type[count];
}
