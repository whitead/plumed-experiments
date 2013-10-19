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

void PREFIX position_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int  i, j, iat;
  int  dir;
  rvec v0, v1;
  rvec pos2;
  real l0, l1, l2;
  int i1, i2, i3;

  for(iat=0;iat<colvar.natoms[i_c];iat++)
    colvar.myder[i_c][iat][0] = colvar.myder[i_c][iat][1] = colvar.myder[i_c][iat][2] = 0.;

  if(colvar.intpar[i_c][0]<3){
    iat = colvar.cvatoms[i_c][0];
    dir = colvar.intpar[i_c][0];
    colvar.ss0[i_c] = 0.;
    for(i=0;i<colvar.natoms[i_c];i++){
      iat = colvar.cvatoms[i_c][i];
      colvar.ss0[i_c] +=  mtd_data->pos[iat][dir]/colvar.natoms[i_c];
      colvar.myder[i_c][i][dir] = 1./colvar.natoms[i_c];
    }

  } else if(colvar.intpar[i_c][0]==3) {
    i1=colvar.intpar[i_c][1];
    i2=colvar.intpar[i_c][2];
    i3=colvar.intpar[i_c][3];
    pos2[0]=pos2[1]=pos2[2]=0.;
    v0[0]=v0[1]=v0[2]=0.;
    v1[0]=v1[1]=v1[2]=0.;
    for(i=0;i<colvar.natoms[i_c];i++){
      iat = colvar.cvatoms[i_c][i];
      for(dir=i1;dir<i2;dir+=i3)pos2[dir]+=mtd_data->pos[iat][dir]/colvar.natoms[i_c];
    }

    for(dir=i1;dir<i2;dir+=i3)v0[dir]=colvar.vecpar[i_c][1][dir]-colvar.vecpar[i_c][0][dir];
    for(dir=i1;dir<i2;dir+=i3)v1[dir]=                 pos2[dir]-colvar.vecpar[i_c][0][dir];
    l0=norm(v0);
    colvar.ss0[i_c]=iprod(v0,v1)/l0;
    for(i=0;i<colvar.natoms[i_c];i++)
      for(dir=i1;dir<i2;dir+=i3)
        colvar.myder[i_c][i][dir]=v0[dir]/colvar.natoms[i_c]/l0;

  } else if(colvar.intpar[i_c][0]==4) {
    i1=colvar.intpar[i_c][1];
    i2=colvar.intpar[i_c][2];
    i3=colvar.intpar[i_c][3];
    pos2[0]=pos2[1]=pos2[2]=0.;
    v0[0]=v0[1]=v0[2]=0.;
    v1[0]=v1[1]=v1[2]=0.;
    for(i=0;i<colvar.natoms[i_c];i++){
      iat = colvar.cvatoms[i_c][i];
      for(dir=i1;dir<i2;dir+=i3)pos2[dir]+=mtd_data->pos[iat][dir]/colvar.natoms[i_c];
    }

    for(dir=i1;dir<i2;dir+=i3)v0[dir]=colvar.vecpar[i_c][1][dir]-colvar.vecpar[i_c][0][dir];
    for(dir=i1;dir<i2;dir+=i3)v1[dir]=                 pos2[dir]-colvar.vecpar[i_c][0][dir];
    l0=norm(v0);
    l1=norm2(v1);
    l2=iprod(v0,v1)/l0;
    colvar.ss0[i_c]=sqrt(l1-l2*l2);
    for(i=0;i<colvar.natoms[i_c];i++)
      for(dir=i1;dir<i2;dir+=i3)
        colvar.myder[i_c][i][dir]=(v1[dir]-l2*v0[dir]/l0)/colvar.natoms[i_c]/colvar.ss0[i_c];

  }
}

// ------------------------------------------------------------------------------------------------

int PREFIX read_position(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw, j;
  double delta = 0.0;
  char string[400];
  char chr[3];
  int help;
  int i1, i2, i3;
  double rtmp;

  help=0;
  chr[0]=chr[1]=chr[2]=' ';

  iw = seek_word(word,"LIST");
  if(iw>=0){   
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR POSITION\n"); help=1;}

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }	

  iw=seek_word(word,"DIR");
  if(iw>=0) { 
    sscanf(word[iw+1],"%s", chr); 
    if(!strcmp(chr,"X") || !strcmp(chr,"x"))
      colvar.intpar[count][0] = 0;
    else if(!strcmp(chr,"Y") || !strcmp(chr,"y"))
      colvar.intpar[count][0] = 1;
    else if(!strcmp(chr,"Z") || !strcmp(chr,"z"))
      colvar.intpar[count][0] = 2;
  }

  iw=seek_word(word,"LINE_POS");
  if(iw>=0) { 
    colvar.intpar[count][0] = 3;
    sscanf(word[++iw],"%s", chr); 
    if     (!strcmp(chr,"X")){
      colvar.intpar[count][1] = 0;
      colvar.intpar[count][2] = 1;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"Y")){
      colvar.intpar[count][1] = 1;
      colvar.intpar[count][2] = 2;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"Z")){
      colvar.intpar[count][1] = 2;
      colvar.intpar[count][2] = 3;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"XY")){
      colvar.intpar[count][1] = 0;
      colvar.intpar[count][2] = 2;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"XZ")){
      colvar.intpar[count][1] = 0;
      colvar.intpar[count][2] = 3;
      colvar.intpar[count][3] = 2;
    }else if(!strcmp(chr,"YZ")){
      colvar.intpar[count][1] = 1;
      colvar.intpar[count][2] = 3;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"XYZ")){
      colvar.intpar[count][1] = 0;
      colvar.intpar[count][2] = 3;
      colvar.intpar[count][3] = 1;
    }else{
      help=1;
    }

    i1=colvar.intpar[count][1];
    i2=colvar.intpar[count][2];
    i3=colvar.intpar[count][3];
    for(i=i1;i<i2;i+=i3){
      sscanf(word[++iw],"%lf", &rtmp); 
      colvar.vecpar[count][0][i]= (real) rtmp; 
    }
    for(i=i1;i<i2;i+=i3){
      sscanf(word[++iw],"%lf", &rtmp); 
      colvar.vecpar[count][1][i]= (real) rtmp; 
    }
  } 

  iw=seek_word(word,"LINE_DIST");
  if(iw>=0) { 
    colvar.intpar[count][0] = 4;
    sscanf(word[++iw],"%s", chr); 
    if     (!strcmp(chr,"X")){
      colvar.intpar[count][1] = 0;
      colvar.intpar[count][2] = 1;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"Y")){
      colvar.intpar[count][1] = 1;
      colvar.intpar[count][2] = 2;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"Z")){
      colvar.intpar[count][1] = 2;
      colvar.intpar[count][2] = 3;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"XY")){
      colvar.intpar[count][1] = 0;
      colvar.intpar[count][2] = 2;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"XZ")){
      colvar.intpar[count][1] = 0;
      colvar.intpar[count][2] = 3;
      colvar.intpar[count][3] = 2;
    }else if(!strcmp(chr,"YZ")){
      colvar.intpar[count][1] = 1;
      colvar.intpar[count][2] = 3;
      colvar.intpar[count][3] = 1;
    }else if(!strcmp(chr,"XYZ")){
      colvar.intpar[count][1] = 0;
      colvar.intpar[count][2] = 3;
      colvar.intpar[count][3] = 1;
    }else{
      help=1;
    }

    i1=colvar.intpar[count][1];
    i2=colvar.intpar[count][2];
    i3=colvar.intpar[count][3];
    for(i=i1;i<i2;i+=i3){
      sscanf(word[++iw],"%lf", &rtmp); 
      colvar.vecpar[count][0][i]= (real) rtmp; 
    }
    for(i=i1;i<i2;i+=i3){
      sscanf(word[++iw],"%lf", &rtmp); 
      colvar.vecpar[count][1][i]= (real) rtmp; 
    }
  } 

  if(chr[0]==' ')help=1;

  if(help){
          fprintf(fplog, "\n-POSIION CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "      POSITION LIST 42 {TYPE} SIGMA 1.0 \n");
          fprintf(fplog, "or in case of groups    \n");
          fprintf(fplog, "      POSITION LIST <g1> {TYPE} SIGMA 1.0 \n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "         6 10    \n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "                 \n");
          fprintf(fplog, " TYPE can be either DIR, LINE_POS or LINE_DIST\n");
          fprintf(fplog, " possible flags for DIR are X, Y and Z\n");
          fprintf(fplog, " possible flags for LINE_POS and LINE_DIST are X, Y, Z, XY, XZ, YZ or XYZ\n");
          plumed_error("PluMeD dead with errors: check log file");
  }

  colvar.type_s[count]   = 32;

  snew(colvar.myder[count], colvar.natoms[count]);

  fprintf(fplog, "%i-%1s POSITION OF %i ATOMS; ", count+1,chr,colvar.natoms[count]);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  fprintf(fplog,"|- SET MEMBERS: ");
  for(i=0;i<colvar.natoms[count];i++){fprintf(fplog," %d ",colvar.cvatoms[count][i]+1);if((i+1)%20==0)fprintf(fplog,"\n               ");}fprintf(fplog,"\n\n");

  return 0;
}
