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

int polycoord_step;

void PREFIX polycoord_newlist(int i_c, struct mtd_data_s *mtd_data)
{
  int firstAtom,secondAtom,i,j,iat;

//  fprintf(stderr,"coord_newslist has been called on step %d.\n",polycoord_step);
  for (j=0; j<colvar.natoms[i_c]; j++) 
  { 
    i = colvar.cvatoms[i_c][j];
    if(j<colvar.list[i_c][0]) (nlist[i_c]).nn[j]=0; 
    (nlist[i_c]).base[j][0]=mtd_data->pos[i][0]; 
    (nlist[i_c]).base[j][1]=mtd_data->pos[i][1]; 
    (nlist[i_c]).base[j][2]=mtd_data->pos[i][2];
  }

  
  rvec rij;
  real mod_rij;
  for(i=0;i<colvar.list[i_c][0];i++) {                                           // sum over CoordNumber(i)
    firstAtom = colvar.cvatoms[i_c][i];
    for(j=colvar.list[i_c][0];j<colvar.natoms[i_c];j++) {
      secondAtom = colvar.cvatoms[i_c][j];
      if(firstAtom == secondAtom) continue;
      if(colvar.cell_pbc[i_c]){
        minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
      } else {
        rij[0] = mtd_data->pos[firstAtom][0]-mtd_data->pos[secondAtom][0];
        rij[1] = mtd_data->pos[firstAtom][1]-mtd_data->pos[secondAtom][1];
        rij[2] = mtd_data->pos[firstAtom][2]-mtd_data->pos[secondAtom][2];
        mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      };
      if (mod_rij<(nlist[i_c]).rskin)
      {
       (nlist[i_c]).ni[i][(nlist[i_c]).nn[i]]=j;(nlist[i_c]).nn[i]++;
       if ((nlist[i_c]).nn[i]>=MAXNN) 
       {
         char buf[1024];
         sprintf(buf,"OOPS! Number of neighbours exceeds MAXNN! Change the source and retry! %d \n",(nlist[i_c]).nn[i]);
         plumed_error(buf); 
       }
      }	  
    }
//    fprintf(stderr," neig:  %d,%d \n",i,(nlist[i_c]).nn[i]);
  }
}

void PREFIX polycoord_checklist(int i_c, struct mtd_data_s *mtd_data)
{
  int i,j; double dr=((nlist[i_c]).rskin-(nlist[i_c]).rcut)*0.5;
  rvec rij;
  real mod_rij;
  for (j=0; j<colvar.natoms[i_c]; j++) 
  { 
    i = colvar.cvatoms[i_c][j];
    if(colvar.cell_pbc[i_c]){
      minimal_image(mtd_data->pos[i], (nlist[i_c]).base[j], &mod_rij, rij);    //UGLY AS HELL, BUT IT'S DLPOLY'S FAULT
    } else {
      rij[0] = mtd_data->pos[i][0]-(nlist[i_c]).base[j][0];
      rij[1] = mtd_data->pos[i][1]-(nlist[i_c]).base[j][1];
      rij[2] = mtd_data->pos[i][2]-(nlist[i_c]).base[j][2];
      mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
    };
    if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { coord_newlist(i_c,mtd_data); break; }
  } 
}

void PREFIX polycoord_restraint_nlist(int i_c, struct mtd_data_s *mtd_data)
{

  int firstAtom, secondAtom, i, j, ix;
  int nn = colvar.nn[i_c], mm = colvar.mm[i_c];

  if(colvar.groups[i_c] == -1) { //only check if we own the list and aren't sharing
    if(polycoord_step==0) coord_newlist(i_c,mtd_data);
    ++polycoord_step; coord_checklist(i_c,mtd_data); 
  }
  
  //this is the possibly shared neighbor list index
  int sn = colvar.groups[i_c] == -1 ? i_c : colvar.groups[i_c];

  rvec rij;
  real mod_rij, rdist, func, dfunc;
  real ncoord, r_0 = colvar.r_0[i_c], d_0 = colvar.d_0[i_c];
  real moment;

  ncoord = 0.;                                                                    
  for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] = 0.;

  for(i=0;i<colvar.list[i_c][0];i++) {                                           // sum over CoordNumber(i)
    firstAtom = colvar.cvatoms[i_c][i];
    for(j=0;j<(nlist[sn]).nn[i];j++) {
//      secondAtom = (nlist[sn]).ni[i][j];
      secondAtom = colvar.cvatoms[i_c][(nlist[sn]).ni[i][j]];
      if(colvar.cell_pbc[i_c]){
        minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
      } else {
        rij[0] = mtd_data->pos[firstAtom][0]-mtd_data->pos[secondAtom][0];
        rij[1] = mtd_data->pos[firstAtom][1]-mtd_data->pos[secondAtom][1];
        rij[2] = mtd_data->pos[firstAtom][2]-mtd_data->pos[secondAtom][2];
        mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      };
      rdist = (mod_rij-d_0)/r_0;
      moment = pow(rdist, colvar.moment[i_c]);
      /* analytic limit of the switching function */
      if(rdist<=0.){
	ncoord+=1;
	dfunc=0;
      }else if(rdist>1){
       dfunc=0.;
      }else{
	func = 2 * rdist * rdist * rdist - 3 * rdist * rdist + 1;
        ncoord += moment * func;
        dfunc = (6 * rdist * rdist - 6 * rdist) / r_0;
	dfunc = (colvar.moment[i_c] * moment * func / rdist + moment * dfunc) / mod_rij;	
      }
      for(ix=0;ix<3;ix++) {
        colvar.myder[i_c][i][ix] += colvar.cn_scale[i_c] * dfunc*rij[ix];
        colvar.myder[i_c][(nlist[sn]).ni[i][j]][ix] += -colvar.cn_scale[i_c] * dfunc*rij[ix];
      }
    }
  }

  colvar.ss0[i_c] = colvar.cn_scale[i_c] * ncoord;

}
void PREFIX polycoord_restraint_no_nlist(int i_c, struct mtd_data_s *mtd_data)
{
  int firstAtom, secondAtom, i, j, ix;
  int nn = colvar.nn[i_c], mm = colvar.mm[i_c];
  rvec rij;
  real mod_rij, rdist, func, dfunc;
  real ncoord, r_0 = colvar.r_0[i_c], d_0 = colvar.d_0[i_c];
  real moment;

  ncoord = 0.;                                                                    
  for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] = 0.;

  for(i=0;i<colvar.list[i_c][0];i++) {                                           // sum over CoordNumber(i)
    firstAtom = colvar.cvatoms[i_c][i];
    for(j=colvar.list[i_c][0];j<colvar.natoms[i_c];j++) {
      if(colvar.logic[i_c] == 1 && (j-colvar.list[i_c][0]) != i) continue;
      secondAtom = colvar.cvatoms[i_c][j];      
      if(firstAtom == secondAtom) continue;
      if(colvar.cell_pbc[i_c]){
        minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
      } else {
        rij[0] = mtd_data->pos[firstAtom][0]-mtd_data->pos[secondAtom][0];
        rij[1] = mtd_data->pos[firstAtom][1]-mtd_data->pos[secondAtom][1];
        rij[2] = mtd_data->pos[firstAtom][2]-mtd_data->pos[secondAtom][2];
        mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      };
      rdist = (mod_rij-d_0)/r_0;
      moment = pow(rdist, colvar.moment[i_c]);
      /* analytic limit of the switching function */
      if(rdist<=0.){
	ncoord+=1;;
	dfunc=0;
      }else if(rdist>1){
       dfunc=0.;
      }else{
	func = 2 * rdist * rdist * rdist - 3 * rdist * rdist + 1;
        ncoord += moment * func;
        dfunc = (6 * rdist * rdist - 6 * rdist) / r_0;
	dfunc = (colvar.moment[i_c] * moment * func / rdist + moment * dfunc) / mod_rij;	
      }
      for(ix=0;ix<3;ix++) {
        colvar.myder[i_c][i][ix] += colvar.cn_scale[i_c] * dfunc*rij[ix];
        colvar.myder[i_c][j][ix] += -colvar.cn_scale[i_c] * dfunc*rij[ix];
      }
    }
  }
  colvar.ss0[i_c] = colvar.cn_scale[i_c] * ncoord;

}

void PREFIX polycoord_restraint(int i_c, struct mtd_data_s *mtd_data) {
  if(logical.nlist[i_c]){
    polycoord_restraint_nlist(i_c,mtd_data);
  }
  else{
    polycoord_restraint_no_nlist(i_c,mtd_data);
  }
}
//--------------------------------------------------------------------------------------------------------------

int PREFIX read_polycoord(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw, iat, j, help, moment;
  double r_0, d_0,r_skin=0.0;
  double delta = 0.0;
  char string[400];

  real threshold, value;
  help=0;
  d_0=0.;
  moment = 0;

  colvar.cell_pbc[count]=1; // default is PBC
  colvar.moment[count]=0; //default
  colvar.b_scale_cn[count] = 0; //default, do not scale
  colvar.cn_scale[count] = 1; //default, do not scale
  colvar.groups[count] = -1;

  iw = seek_word(word,"LIST");
  if(iw>=0){   
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR POLYCOORD\n"); help=1;}
  iw=seek_word(word,"PAIR");
  if(iw>=0) colvar.logic[count] = 1;
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.mm[count]);} else { fprintf(fplog,"|- NEEDED MM KEYWORD FOR POLYCOORD\n"); help=1;}
  iw=seek_word(word,"MOMENT");
  if(iw>=0) { sscanf(word[iw+1],"%i", &moment);}
  iw=seek_word(word,"SHARE_NLIST");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.groups[count]);}
  iw=seek_word(word,"SCALE");
  if(iw>=0) { colvar.b_scale_cn[count] = 1;}
  iw=seek_word(word,"R_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &r_0); } else { fprintf(fplog,"|- NEEDED R_0 KEYWORD FOR POLYCOORD\n"); help=1;}
  iw=seek_word(word,"D_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &d_0); }
  iw=seek_word(word,"R_SKIN");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &r_skin);} 
  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}

  if(help){
          fprintf(fplog, "\n-POLYCOORD CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "POLYCOORD LIST <g1> <g2>  NN 6 MM 12 R_0 0.75 D_0 3. SIGMA 1.0 \n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "         5 1 6    \n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "       \n");
          fprintf(fplog, "         g2->    \n");
          fprintf(fplog, "         LOOP 6921 21786 5 \n");
          fprintf(fplog, "         g2<-    \n");
          fprintf(fplog, "       \n");
          plumed_error("PluMed dead with errors: check log file");
  }

  colvar.r_0[count]      = (real) r_0;
  colvar.d_0[count]      = (real) d_0;
  colvar.type_s[count]   = 3;
  colvar.moment[count]   = moment;

  fprintf(fplog, "%1i-POLYCOORDINATION NUMBER OF (1st SET: %i ATOMS) WRT (2nd SET: %i ATOMS); ", count+1, colvar.list[count][0], colvar.list[count][1]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON");
  else                       fprintf(fplog, " PBC OFF");
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  fprintf(fplog, "|--FUNCTIONAL FORM: r^(moment) * (1-((d_ij-d_0)/r_0)^n) / (1-((d_ij-d_0)/r_0)^m) \n");
  fprintf(fplog, "|--PARAMETERS: r_0= %f d_0= %f moment= %i\n", colvar.r_0[count], colvar.d_0[count], colvar.moment[count]);
    

  if(colvar.logic[count]) fprintf(fplog,"|--PAIR POLYCOORDINATION ACTIVE \n"); 
  if(colvar.logic[count] == 1 && colvar.list[count][0] != colvar.list[count][1])
    plumed_error("If PAIR is active <g1> anf <g2> must have the same number of atoms");


  snew(colvar.myder[count], colvar.natoms[count]);

  //allocate arrays for internal Verlet lists machinery
  if(r_skin>0.0)  
  {  
    if(colvar.groups[count] > -1) plumed_error("CAN ONLY USE EITHER R_SKIN OR SHARE_NLIST");
     if (colvar.logic[count]) plumed_error("PAIR OPTION IS NOT WORKING WITH VERLET LIST");
     fprintf(fplog, "|--VERLET LIST ACTIVE: R_CUT %lf R_SKIN %lf \n",r_0+d_0,r_skin); 
     logical.nlist[count]=1;
     (nlist[count]).rskin=r_skin;  
     (nlist[count]).rcut=r_0 + d_0; 
     (nlist[count]).base=float_2d_array_alloc(colvar.natoms[count],3);
     (nlist[count]).nn=(int *) malloc(colvar.list[count][0]*sizeof(int));
     (nlist[count]).ni=int_2d_array_alloc(colvar.list[count][0],MAXNN);
     polycoord_step=0;
  } else if(colvar.groups[count] > -1) {
    colvar.groups[count] -= 1; //convert index from input (1->) to code (0->)
    fprintf(fplog,"|--WILL SHARE NLIST WITH CV %d\n", colvar.groups[count] + 1);    
    logical.nlist[count] = 1;
  } else
  {fprintf(fplog,"|--VERLET LIST NOT ACTIVE\n");}
  
  if(colvar.b_scale_cn[count]) {
    colvar.cn_scale[count] = 1. / colvar.list[count][0];
    fprintf(fplog, "|--WILL SCALE POLYCOORDINATION NUMBER BY ATOM NUMBER (%d)--|\n", colvar.list[count][0] );
  }

  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  return colvar.natoms[count]; 
}

