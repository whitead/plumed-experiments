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


void PREFIX poly_restraint(int i_c, struct mtd_data_s *mtd_data)
{

int i,ix;
double shift, coeff, sterm, expo, vcv, factor;
int jcv;
int iatG, iatT;

// fprintf(mtd_data->fplog,"count is %d\n",i_c);
// for(int i=0;i<i_c; i++) fprintf(mtd_data->fplog,"%d, num=%f\n", i, colvar.ss0[i]);

i=0;
colvar.ss0[i_c]=0;
if(1) while(colvar.intpar[i_c][i]>=0){
   jcv = colvar.intpar[i_c][i];
   shift = colvar.vecpar[i_c][i][0];
   coeff=  colvar.vecpar[i_c][i][1];
   expo =  colvar.vecpar[i_c][i][2];
   sterm = coeff*pow(colvar.ss0[jcv]+shift,expo);
    colvar.ss0[i_c] += sterm;
   /* fprintf(mtd_data->fplog, 
     "| poly TERM %d (cv=%d, sh=%f, coef=%f, exp=%f) CONTRIB %f\n", 
     i, jcv, shift, coeff, expo, sterm); */
   i++;
}

//fprintf(mtd_data->fplog,"natoms=%d\n",colvar.natoms[i_c]);
for(i=0;i<colvar.natoms[i_c];i++){ 
 for(ix=0;ix<3;ix++){ 
  colvar.myder[i_c][i][ix]=0; } }


// Calculates derivative
iatG=0;
i=0;
while(colvar.intpar[i_c][i]>=0){
  jcv=colvar.intpar[i_c][i];
  shift = colvar.vecpar[i_c][i][0];
  coeff=  colvar.vecpar[i_c][i][1];
  expo =  colvar.vecpar[i_c][i][2];
  if(fabs(expo-0.00)<1e-7) {vcv=0.0;} 
   else if(fabs(expo-1.0)<1e-7) { vcv = 1.0;} 
   else { vcv = pow(colvar.ss0[jcv]+shift,expo-1);} 
  factor=expo*coeff*vcv;
  for(iatT=0;iatT<colvar.natoms[jcv];iatT++){
    for(ix=0;ix<3;ix++) colvar.myder[i_c][iatG][ix] += factor*colvar.myder[jcv][iatT][ix];
     // for debug
     /* fprintf(mtd_data->fplog, 
       "| poly(%d) term=%d (iatG=%d [%5.3f,%5.3f,%5.3f], fctr=%f, jcv=%d iatT=%d [%5.3f,%5.3f,%5.3f] \n", 
        i, i_c, iatG, colvar.myder[i_c][iatG][0], colvar.myder[i_c][iatG][1], 
                 colvar.myder[i_c][iatG][2], 
        factor, jcv, iatT, 
        colvar.myder[jcv][iatT][0], colvar.myder[jcv][iatT][1], 
        colvar.myder[jcv][iatT][2]
    ); */
    iatG++;
  }
  i++; 
}


return;
}



int PREFIX read_poly(char **word, int count, t_plumed_input *input,int *iline,FILE *fplog)
{

int iw, i, j, k;
double delta=0.0;
int nterms;
int help;
help=0;
double shift, coeff,expo ;
int jcv;

// fprintf(fplog,"count is %d\n",count);
// for(int i=0;i<count; i++) fprintf(fplog,"%d, num=%d\n", i, colvar.natoms[i]);

iw = seek_word(word,"TERMS");
if(iw>=0) { sscanf(word[iw+1],"%d", &nterms); } 
 else { fprintf(fplog,"|- NEEDED TERMS KEYWORD FOR POLY\n"); help=1;}
if(nterms>20) {
  fprintf(fplog,"|- TO USE MORE THAN 20 TERMS, MODIFY ALLOCATION OF vecpar in metadyn.h!\n");
  plumed_error("PluMeD dead with errors: check log file");
  exit(1); // just in case...
}

iw=seek_word(word,"SIGMA");
if(iw>=0){ 
  sscanf(word[iw+1],"%lf", &delta);
  colvar.delta_r[count]  = (real) delta; 
}

 if(help){
  fprintf(fplog, "\n-POLY CV: WRONG SYNTAX\n");
  fprintf(fplog, "e.g.:     \n");
  fprintf(fplog, "      POLY TERMS 3 {SIGMA 1.0}\n");
  fprintf(fplog, "      CV 1 {SHIFT 0.0} {COEFF 1.0} {EXP 1.0}\n");
  fprintf(fplog, "      CV 2 {SHIFT 0.0} {COEFF 1.0} {EXP 1.0}\n");
  plumed_error("PluMeD dead with errors: check log file");
 }


fprintf(fplog, "\n%li-POLY (POLYNOMIAL FUNCTION OF CVs): COMBINING %d TERMS, ",count+1,nterms);
if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
 else fprintf(fplog,"\n");

k=0;
colvar.natoms[count]=0;
for(i=0; i<nterms; i++){
    fprintf(fplog, "|- TERM %d ",i+1);
    shift = 0.0;
    coeff = 1.0;
    expo  = 1.0;
    jcv = 0;
    (*iline)++;
    for(j=0;j<(input->nwords[*iline])-1;j++) {
     if(seek_word(&(input->words[*iline][j]),"CV")>=0){
      sscanf(input->words[*iline][j+1],"%d",&jcv);}
     if(seek_word(&(input->words[*iline][j]),"SHIFT")>=0){
       sscanf(input->words[*iline][j+1],"%lf",&shift); }
     if(seek_word(&(input->words[*iline][j]),"COEFF")>=0){ 
       sscanf(input->words[*iline][j+1],"%lf",&coeff); }
     if(seek_word(&(input->words[*iline][j]),"EXP")>=0){ 
       sscanf(input->words[*iline][j+1],"%lf",&expo); }
    }
    if(jcv==0){
     fprintf(fplog,"\n\nCV NEEDS TO BE SPECIFIED IN POLY!\n");
     plumed_error("PluMeD dead with errors: check log file");
     exit(1); // just to be sure!
    } else if (jcv>count) {
     fprintf(fplog,"\n\nCV %d DOES NOT APPEAR TO HAVE BEEN DEFINED YET.\n",jcv+1);
     fprintf(fplog,"ONLY PREVIOUSLY DEFINED CVs CAN BE COMBINED BY POLY!\n");
     plumed_error("PluMeD dead with errors: check log file");
     exit(1); // just to be sure!
    } else {
     jcv--;
     colvar.intpar[count][i] = jcv;
     logical.always[ jcv ]=1;     
     colvar.vecpar[count][i][0] = shift;
     colvar.vecpar[count][i][1] = coeff;
     colvar.vecpar[count][i][2] = expo;
     colvar.natoms[count]+=colvar.natoms[jcv]; 
     // fprintf(fplog,"jcv=%d, natoms=%d\n", jcv, colvar.natoms[jcv]);
     srenew(colvar.cvatoms[count],colvar.natoms[count]);
     // Copy atom indexes from jcv into count 
     for(j=0;j<colvar.natoms[jcv];j++){
         colvar.cvatoms[count][k] = colvar.cvatoms[jcv][j]; 
         k++; }
    }
    fprintf(fplog," CV %d SHIFT %lf, COEFF %lf, EXP %lf (NATOMS=%d)\n", jcv, shift,coeff,expo, 
    colvar.natoms[count]);
}

 // a negative value in intpar means no more terms
 colvar.intpar[count][nterms] = -1;


// FOR DEBUG
if(0){ 
 i=0;
 while(colvar.intpar[count][i]>=0){
  fprintf(fplog,"|TERM %d: CV %d SHIFT %lf, COEFF %lf, EXP %lf\n", 
   i, colvar.intpar[count][i],colvar.vecpar[count][i][0], 
      colvar.vecpar[count][i][1],colvar.vecpar[count][i][2]); 
   i++;
 }
 fprintf(fplog,"|%d ATOMS IN THIS CV:\n", colvar.natoms[count]);
 for(j=0; j<colvar.natoms[count]; j++)
  fprintf(fplog, "|-(ipoly=%d), (atid=%d)\n",j,colvar.cvatoms[count][j]);
} 

  colvar.type_s[count]   = 50;
  snew(colvar.myder[count], colvar.natoms[count]);
  return colvar.natoms[count];

}


