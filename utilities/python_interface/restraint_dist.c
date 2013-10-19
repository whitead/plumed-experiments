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
void PREFIX dist_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i, iat;
  real mod_rij, mass1, mass2;
  rvec rij, sum1, sum2;
	
  mass1 = mass2 = 0.;
  sum1[0] = sum1[1] = sum1[2] = 0.;
  sum2[0] = sum2[1] = sum2[2] = 0.;


  // point to axis
  if(colvar.intpar[i_c][3]==1){
     pt_from_axis_restraint(i_c, mtd_data);  
     return;
  } else if(colvar.intpar[i_c][3]==2){
     diffdist_restraint(i_c, mtd_data);  
     return;
  } else if(colvar.intpar[i_c][7]==1){
     proj_on_axis_restraint(i_c, mtd_data);
     return;
  };
  
  for(i=0;i<colvar.list[i_c][0];i++) {
    iat = colvar.cvatoms[i_c][i];
    sum1[0] += mtd_data->mass[iat]*mtd_data->pos[iat][0];
    sum1[1] += mtd_data->mass[iat]*mtd_data->pos[iat][1];
    sum1[2] += mtd_data->mass[iat]*mtd_data->pos[iat][2];
    mass1 += mtd_data->mass[iat];
  }
  for(i=colvar.list[i_c][0];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    sum2[0] += mtd_data->mass[iat]*mtd_data->pos[iat][0];
    sum2[1] += mtd_data->mass[iat]*mtd_data->pos[iat][1];
    sum2[2] += mtd_data->mass[iat]*mtd_data->pos[iat][2];
    mass2 += mtd_data->mass[iat];
  }

  sum1[0] /= mass1; sum1[1] /= mass1; sum1[2] /= mass1;
  sum2[0] /= mass2; sum2[1] /= mass2; sum2[2] /= mass2;

  // Project center of masses, if required
  if(colvar.intpar[i_c][0])
    sum2[0]=sum1[0]=0.0;
  if(colvar.intpar[i_c][1])
    sum2[1]=sum1[1]=0.0;
  if(colvar.intpar[i_c][2])
    sum2[2]=sum1[2]=0.0;


  if(colvar.cell_pbc[i_c]){
    minimal_image(sum1, sum2, &mod_rij, rij);
  } else {
    rij[0] = sum1[0]-sum2[0];
    rij[1] = sum1[1]-sum2[1];
    rij[2] = sum1[2]-sum2[2];
    mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
  };

  for(i=0;i<colvar.list[i_c][0];i++) { 
    iat = colvar.cvatoms[i_c][i];
    colvar.myder[i_c][i][0] =  (mtd_data->mass[iat]*rij[0])/(mass1*mod_rij);
    colvar.myder[i_c][i][1] =  (mtd_data->mass[iat]*rij[1])/(mass1*mod_rij);
    colvar.myder[i_c][i][2] =  (mtd_data->mass[iat]*rij[2])/(mass1*mod_rij);
  }
  for(i=colvar.list[i_c][0];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    colvar.myder[i_c][i][0] = -(mtd_data->mass[iat]*rij[0])/(mass2*mod_rij);
    colvar.myder[i_c][i][1] = -(mtd_data->mass[iat]*rij[1])/(mass2*mod_rij);
    colvar.myder[i_c][i][2] = -(mtd_data->mass[iat]*rij[2])/(mass2*mod_rij);
  }

  colvar.ss0[i_c] = mod_rij;

}

// ------------------------------------------------------------------------------------------------

int PREFIX read_dist(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw, iat,j;
  double delta = 0.0;
  char string[400];
  char *name_ow1=NULL,*name_ow2=NULL,*name_ow3=NULL;
  int help;
  help=0;

  colvar.cell_pbc[count]=1; // default is PBC

  // default is not to project
  colvar.intpar[count][0] = 0;
  colvar.intpar[count][1] = 0;
  colvar.intpar[count][2] = 0;
  // point to axis is by default off
  colvar.intpar[count][3] = 0;
  // optimized external weights
  colvar.intpar[count][4] = 0;
  colvar.intpar[count][5] = 0;
  colvar.intpar[count][6] = 0;
  // proj to axis is by default off
  colvar.intpar[count][7] = 0;

  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}

  iw = seek_word(word,"LIST");
  if(iw>=0){   
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR DISTANCE\n"); help=1;}

  /* Project on a given axis. (TONI) I am putting this before SIGMA
     for consisency with POSITION  */
  iw=seek_word(word,"DIR");
  if(iw>=0) { 
    char chr[3];
    int dropx,dropy,dropz;	/* discard x,y,z component? */
    dropx=dropy=dropz=1;	/* yes on all axis, except those requested */
    sscanf(word[iw+1],"%s", chr); 
    if(strchr(chr,'X') || strchr(chr,'x')) {
      dropx=0;			/* keep x */
    } 
    if(strchr(chr,'Y') || strchr(chr,'y')) {
      dropy=0;
    }
    if(strchr(chr,'Z') || strchr(chr,'z')) {
      dropz=0;
    }
    if(dropx && dropy && dropz) {
      fprintf(fplog,"|- INVALID PROJECTION SPECIFIED\n"); help=1;
    } 
    colvar.intpar[count][0] = dropx;
    colvar.intpar[count][1] = dropy;
    colvar.intpar[count][2] = dropz;
  }

  /* point from axis : in this case take another point or group: the one that appear in list define the axis, this defines the point   */
  iw = seek_word(word,"POINT_FROM_AXIS");
  if(iw>=0){   
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][2]=j;
             /* switch on point from axis option */
             colvar.intpar[count][3] = 1;
  }
  /* proj on axis : in this case take another point or group: the one that appear in list define the axis, this defines the point   */
  iw = seek_word(word,"PROJ_ON_AXIS");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][2]=j;
             /* switch on proj on axis option */
             colvar.intpar[count][7] = 1;
  } 
  // optimized external weights only for this variable 
  iw = seek_word(word,"OW1");
  if(iw>=0){   
      colvar.intpar[count][4]=1; 
      name_ow1=word[iw+1]; 
  }
  iw = seek_word(word,"OW2");
  if(iw>=0){   
      colvar.intpar[count][5]=1; 
      name_ow2=word[iw+1]; 
  }
  iw = seek_word(word,"OW3");
  if(iw>=0){   
      colvar.intpar[count][6]=1; 
      name_ow3=word[iw+1]; 
  }
  // other two groups if you need difference of distances  
  iw = seek_word(word,"DIFFDIST");
  if(iw>=0){   
        j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
        colvar.natoms[count]+=j;
        colvar.list[count][2]=j;
        j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
        colvar.natoms[count]+=j;
        colvar.list[count][3]=j;
        /* switch on  diffdist option */
        colvar.intpar[count][3] = 2;
  } 
 
 
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);  
             colvar.delta_r[count]  = (real) delta; }

 // else if (logical.do_hills) { fprintf(fplog,"|- NEEDED SIGMA KEYWORD FOR DISTANCE\n"); help=1;}
  if(help){
          fprintf(fplog, "\n-DISTANCE CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "      DISTANCE LIST 15 30 DIR {projection} { POINT_FROM_AXIS/PROJ_ON_AXIS <g3>  }  {DIFFDIST  <g3> <g4> }SIGMA 1.0 \n");
          fprintf(fplog, "  \n");
          fprintf(fplog, "or in case of groups    \n");
          fprintf(fplog, "      DISTANCE LIST <g1> <g2> SIGMA 1.0 \n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "         6 10    \n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "                 \n");
          fprintf(fplog, "         g2->    \n");
          fprintf(fplog, "         8 15 21 \n");
          fprintf(fplog, "         g2<-    \n");
          fprintf(fplog, "                 \n");
          fprintf(fplog, " {projection} can be X, Y, Z, XY, XZ or YZ");
          plumed_error("PluMeD dead with errors: check log file");
  }

  colvar.type_s[count]   = 1;

  fprintf(fplog, "\n%1i-DISTANCE: (1st SET: %i ATOMS), (2nd SET: %i ATOMS); ", count+1, colvar.list[count][0], colvar.list[count][1]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON");
  else                       fprintf(fplog, " PBC OFF");
  if (logical.do_hills){
 	if (colvar.delta_r[count]>0){
        	 fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
        }
  }
  else fprintf(fplog,"\n");

  iat=0;
  fprintf(fplog, "|- DISCARDING DISTANCE COMPONENTS (XYZ): %d%d%d\n",
	  colvar.intpar[count][0],colvar.intpar[count][1],colvar.intpar[count][2]);
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");
  if(colvar.intpar[count][3]==1){
    fprintf(fplog,"|- 3nd SET MEMBERS (ONLY FOR POINT_FROM_AXIS: this defines the point, the others above define the two points of the axis): ");
    for(i=0;i<colvar.list[count][2];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");
  }
  if(colvar.intpar[count][7]==1){
    fprintf(fplog,"|- 3nd SET MEMBERS (ONLY FOR PROJ_ON_AXIS: this defines the point, the others above define the two points of the axis): ");
    for(i=0;i<colvar.list[count][2];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");
  }
  if(colvar.intpar[count][3]==2){
    fprintf(fplog,"|- 3nd SET MEMBERS (ONLY FOR DIFFDIST: this defines the other two points): ");
    for(i=0;i<colvar.list[count][2];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
    fprintf(fplog,"|- 3nd SET MEMBERS (ONLY FOR DIFFDIST: this defines the other two points): ");
    for(i=0;i<colvar.list[count][3];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");
  }


  /* import the optimized weights */
  if(colvar.intpar[count][4]){
     // pass the address to the first  address of region
     import_ow(&(colvar.ow_weight[count*3]),name_ow1,0,colvar.list[count][0],colvar.cvatoms[count]); 
  } 
  if(colvar.intpar[count][5]){
     import_ow(&(colvar.ow_weight[count*3+1]),name_ow2,colvar.list[count][0],colvar.list[count][0]+colvar.list[count][1],colvar.cvatoms[count]); 
  } 
  if(colvar.intpar[count][6]){
     import_ow(&(colvar.ow_weight[count*3+2]),name_ow3,colvar.list[count][1]+colvar.list[count][0],colvar.list[count][0]+colvar.list[count][1]+colvar.list[count][2],colvar.cvatoms[count]); 
  } 

  snew(colvar.myder[count], colvar.natoms[count]);

  return colvar.natoms[count];
}
void PREFIX pt_from_axis_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i,j, iat;
  real mod_ba, mod_ha, mass_a, mass_b, mass_h, dot, mod_ptx, mass;
  rvec r_ba,r_ha, sum_a, sum_b, sum_h;
  real tau_t;

  mass_a = mass_b = mass_h = 0.;
  sum_a[0] = sum_a[1] = sum_a[2] = 0.;
  sum_b[0] = sum_b[1] = sum_b[2] = 0.;
  sum_h[0] = sum_h[1] = sum_h[2] = 0.;

  j=0;
  // added optimized weights for axis fitting
  for(i=0;i<colvar.list[i_c][0];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][4]){mass=colvar.ow_weight[i_c*3][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    sum_a[0] += mass*mtd_data->pos[iat][0];
    sum_a[1] += mass*mtd_data->pos[iat][1];
    sum_a[2] += mass*mtd_data->pos[iat][2];
    if(colvar.intpar[i_c][4]){mass_a +=1;}
    else {mass_a += mass;}
  }
  j=0;
  for(i=colvar.list[i_c][0];i<colvar.list[i_c][0]+colvar.list[i_c][1];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][5]){mass=colvar.ow_weight[i_c*3+1][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    sum_b[0] += mass*mtd_data->pos[iat][0];
    sum_b[1] += mass*mtd_data->pos[iat][1];
    sum_b[2] += mass*mtd_data->pos[iat][2];
    if(colvar.intpar[i_c][5]){mass_b +=1;}
    else {mass_b += mass;}
  }
  j=0;
  for(i=colvar.list[i_c][0]+colvar.list[i_c][1];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][6]){mass=colvar.ow_weight[i_c*3+2][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    sum_h[0] += mass*mtd_data->pos[iat][0];
    sum_h[1] += mass*mtd_data->pos[iat][1];
    sum_h[2] += mass*mtd_data->pos[iat][2];
    if(colvar.intpar[i_c][6]){mass_h +=1;}
    else {mass_h += mass;}
  }


  sum_a[0] /= mass_a; sum_a[1] /= mass_a; sum_a[2] /= mass_a;
  sum_b[0] /= mass_b; sum_b[1] /= mass_b; sum_b[2] /= mass_b;
  sum_h[0] /= mass_h; sum_h[1] /= mass_h; sum_h[2] /= mass_h;

  // 

  if(colvar.cell_pbc[i_c]){
    minimal_image(sum_b, sum_a, &mod_ba, r_ba);
  } else {
    r_ba[0] = sum_b[0]-sum_a[0];
    r_ba[1] = sum_b[1]-sum_a[1];
    r_ba[2] = sum_b[2]-sum_a[2];
    mod_ba  = sqrt(r_ba[0]*r_ba[0]+r_ba[1]*r_ba[1]+r_ba[2]*r_ba[2]);
  };

  if(colvar.cell_pbc[i_c]){
    minimal_image(sum_h, sum_a, &mod_ha, r_ha);
  } else {
    r_ha[0] = sum_h[0]-sum_a[0];
    r_ha[1] = sum_h[1]-sum_a[1];
    r_ha[2] = sum_h[2]-sum_a[2];
    mod_ha  = sqrt(r_ha[0]*r_ha[0]+r_ha[1]*r_ha[1]+r_ha[2]*r_ha[2]);
  };

  // dot
  dot=r_ha[0]*r_ba[0]+r_ha[1]*r_ba[1]+r_ha[2]*r_ba[2];

  // point to axis dist
  mod_ptx=pow(mod_ha*mod_ha-(dot*dot)/(mod_ba*mod_ba),0.5);  

  // tau_t
  tau_t=dot/(mod_ba*mod_ba)  ;

  // put everything in globally visible structures

  j=0;
  for(i=0;i<colvar.list[i_c][0];i++) { 
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][4]){mass=colvar.ow_weight[i_c*3][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    colvar.myder[i_c][i][0] = mass/mass_a*(-r_ha[0]+(r_ba[0]+r_ha[0])*tau_t-r_ba[0]*tau_t*tau_t)/mod_ptx ; 
    colvar.myder[i_c][i][1] = mass/mass_a*(-r_ha[1]+(r_ba[1]+r_ha[1])*tau_t-r_ba[1]*tau_t*tau_t)/mod_ptx ; 
    colvar.myder[i_c][i][2] = mass/mass_a*(-r_ha[2]+(r_ba[2]+r_ha[2])*tau_t-r_ba[2]*tau_t*tau_t)/mod_ptx ;
  }
  j=0;
  for(i=colvar.list[i_c][0];i<colvar.list[i_c][0]+colvar.list[i_c][1];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][5]){mass=colvar.ow_weight[i_c*3+1][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    colvar.myder[i_c][i][0] = mass/mass_b*(-r_ha[0]*tau_t+r_ba[0]*tau_t*tau_t)/mod_ptx ; 
    colvar.myder[i_c][i][1] = mass/mass_b*(-r_ha[1]*tau_t+r_ba[1]*tau_t*tau_t)/mod_ptx ; 
    colvar.myder[i_c][i][2] = mass/mass_b*(-r_ha[2]*tau_t+r_ba[2]*tau_t*tau_t)/mod_ptx ; 
  }
  j=0;
  for(i=colvar.list[i_c][0]+colvar.list[i_c][1];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][6]){mass=colvar.ow_weight[i_c*3+2][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    colvar.myder[i_c][i][0] = mass/mass_h*(r_ha[0]-r_ba[0]*tau_t)/mod_ptx ;
    colvar.myder[i_c][i][1] = mass/mass_h*(r_ha[1]-r_ba[1]*tau_t)/mod_ptx ;
    colvar.myder[i_c][i][2] = mass/mass_h*(r_ha[2]-r_ba[2]*tau_t)/mod_ptx ;
  }

  colvar.ss0[i_c] = mod_ptx;

  if( !(colvar.it%colvar.nt_print) ) {
     fprintf(mtd_data->fplog,"|-POINT_FROM_AXIS: A %8.3f %8.3f %8.3f  ",sum_a[0],sum_a[1],sum_a[2]);
     fprintf(mtd_data->fplog," B %8.3f %8.3f %8.3f  ",sum_b[0],sum_b[1],sum_b[2]);
     fprintf(mtd_data->fplog," T %8.3f %8.3f %8.3f  ",sum_b[0]+tau_t*(sum_b[0]-sum_a[0]),sum_b[1]+tau_t*(sum_b[1]-sum_a[1]),sum_b[2]+tau_t*(sum_b[2]-sum_a[2]));
     fprintf(mtd_data->fplog," H %8.3f %8.3f %8.3f  \n",sum_h[0],sum_h[1],sum_h[2]);
  }



}
void PREFIX import_ow(real **weight ,char* filename , int start, int end , int *indexes )
{
    FILE *fp;
    char string[100],*str; 
    int id,size,i,j;
    int *ind_tmp ;
    int *idone; 
    real *weight_tmp;
    real ww;
    fprintf(mtd_data.fplog,"|-IMPORTING OPTIMIZED MASS WEIGHT %s\n",filename);
    fp=fopen(filename,"r");
    if (fp == NULL) {
       char buf[1024];
       sprintf(buf,"UNSUCCESSFULL OPENING FOR FILE %s",filename);
       plumed_error(buf);  
    }
    j=0;
    while(1){
        str=fgets(string,100,fp);
        if(str==NULL)break;
        if (feof(fp))break;
        j++;
    }
    fclose(fp);
    snew(ind_tmp,j);snew(weight_tmp,j);

    size=0;
    fp=fopen(filename,"r");
    while(1){
        str=fgets(string,100,fp);
        if(str==NULL)break;
        if (feof(fp))break;
        char *result = NULL;
        result = strtok( string, " \t" );   
        id=atoi(result);
        result = strtok( NULL, " \t" );   
        ww=atof(result);
        weight_tmp[size]=ww; ind_tmp[size]=id-1;
        //fprintf(mtd_data.fplog,"|-IND %d WEIGHT %f \n",ind_tmp[size]+1,weight_tmp[size]);
        size++; 
    }
    // check consistency
    // number first 
    if(size!=(end-start)) {
       char buf[1024];
       sprintf(buf,"Expected %d mass weights. Found %d. Only one weight per atom is allowed!",end-start,size);
       plumed_error(buf);  
    }
    snew(*weight,size);
    //  reorder and check 
    snew(idone,size);for(i=0;i<size;i++)idone[i]=0; 
    for(i=0;i<size;i++){
          if(!idone[i]){ // only if it is not done
                for(j=start;j<end;j++){
                     if(indexes[j]==ind_tmp[i]){ 
                        //(*weight)[j]=weight_tmp[i];
                        (*weight)[i]=weight_tmp[i];
                        idone[i]=1;
                        fprintf(mtd_data.fplog,"|-FOUND_MATCH FOR ATOM %d : WW %f \n",indexes[j]+1,(*weight)[i]);
                        break;
                     } 
                } 
          }else{ plumed_error("The same atom is occurring twice in the weights"); }
    }
  sfree(ind_tmp); 
  sfree(weight_tmp); 
  sfree(idone);
  return;
} 
void PREFIX diffdist_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i,j, iat;
  real mod_ba, mod_gh, mass_a, mass_b, mass_h, mass_g, mod_diff,mass;
  rvec r_ba,r_gh, sum_a, sum_b, sum_h, sum_g;
  rvec dd_drh,dd_dra,dd_drb,dd_drg; 

  mass_a = mass_b = mass_h = mass_g= 0.;
  sum_a[0] = sum_a[1] = sum_a[2] = 0.;
  sum_b[0] = sum_b[1] = sum_b[2] = 0.;
  sum_h[0] = sum_h[1] = sum_h[2] = 0.;
  sum_g[0] = sum_g[1] = sum_g[2] = 0.;

  for(i=0;i<colvar.list[i_c][0];i++) {
    iat = colvar.cvatoms[i_c][i];
    mass=mtd_data->mass[iat];   
    sum_a[0] += mass*mtd_data->pos[iat][0];
    sum_a[1] += mass*mtd_data->pos[iat][1];
    sum_a[2] += mass*mtd_data->pos[iat][2];
    mass_a += mass;
  }
  for(i=colvar.list[i_c][0];i<colvar.list[i_c][0]+colvar.list[i_c][1];i++) {
    iat = colvar.cvatoms[i_c][i];
    mass=mtd_data->mass[iat];   
    sum_b[0] += mass*mtd_data->pos[iat][0];
    sum_b[1] += mass*mtd_data->pos[iat][1];
    sum_b[2] += mass*mtd_data->pos[iat][2];
    mass_b += mass;
  }
  for(i=colvar.list[i_c][0]+colvar.list[i_c][1];i<colvar.list[i_c][0]+colvar.list[i_c][1]+colvar.list[i_c][2];i++) {
    iat = colvar.cvatoms[i_c][i];
    mass=mtd_data->mass[iat];   
    sum_h[0] += mass*mtd_data->pos[iat][0];
    sum_h[1] += mass*mtd_data->pos[iat][1];
    sum_h[2] += mass*mtd_data->pos[iat][2];
    mass_h += mass;
  }
  for(i=colvar.list[i_c][0]+colvar.list[i_c][1]+colvar.list[i_c][2];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    mass=mtd_data->mass[iat];   
    sum_g[0] += mass*mtd_data->pos[iat][0];
    sum_g[1] += mass*mtd_data->pos[iat][1];
    sum_g[2] += mass*mtd_data->pos[iat][2];
    mass_g += mass;
  }


  sum_a[0] /= mass_a; sum_a[1] /= mass_a; sum_a[2] /= mass_a;
  sum_b[0] /= mass_b; sum_b[1] /= mass_b; sum_b[2] /= mass_b;
  sum_h[0] /= mass_h; sum_h[1] /= mass_h; sum_h[2] /= mass_h;
  sum_g[0] /= mass_g; sum_g[1] /= mass_g; sum_g[2] /= mass_g;

  // 

  if(colvar.cell_pbc[i_c]){
    minimal_image(sum_b, sum_a, &mod_ba, r_ba);
    minimal_image(sum_g, sum_h, &mod_gh, r_gh);
  } else {
    r_ba[0] = sum_b[0]-sum_a[0];
    r_ba[1] = sum_b[1]-sum_a[1];
    r_ba[2] = sum_b[2]-sum_a[2];
    mod_ba  = sqrt(r_ba[0]*r_ba[0]+r_ba[1]*r_ba[1]+r_ba[2]*r_ba[2]);
    r_gh[0] = sum_g[0]-sum_h[0];
    r_gh[1] = sum_g[1]-sum_h[1];
    r_gh[2] = sum_g[2]-sum_h[2];
    mod_gh  = sqrt(r_gh[0]*r_gh[0]+r_gh[1]*r_gh[1]+r_gh[2]*r_gh[2]);
  };

  // point to axis dist

  mod_diff=mod_ba-mod_gh;  

  // derivative respect to point

  dd_dra[0]= (sum_a[0]-sum_b[0])/mod_ba;
  dd_dra[1]= (sum_a[1]-sum_b[1])/mod_ba;
  dd_dra[2]= (sum_a[2]-sum_b[2])/mod_ba;

  dd_drb[0]= (sum_b[0]-sum_a[0])/mod_ba; 
  dd_drb[1]= (sum_b[1]-sum_a[1])/mod_ba; 
  dd_drb[2]= (sum_b[2]-sum_a[2])/mod_ba; 

  dd_drh[0]= (sum_g[0]-sum_h[0])/mod_gh;
  dd_drh[1]= (sum_g[1]-sum_h[1])/mod_gh;
  dd_drh[2]= (sum_g[2]-sum_h[2])/mod_gh;

  dd_drg[0]= (sum_h[0]-sum_g[0])/mod_gh; 
  dd_drg[1]= (sum_h[1]-sum_g[1])/mod_gh; 
  dd_drg[2]= (sum_h[2]-sum_g[2])/mod_gh; 
 
  // put everything in globally visible structures

  for(i=0;i<colvar.list[i_c][0];i++) { 
    iat = colvar.cvatoms[i_c][i];
    mass=mtd_data->mass[iat];   
    colvar.myder[i_c][i][0] =  (mass*dd_dra[0])/(mass_a);
    colvar.myder[i_c][i][1] =  (mass*dd_dra[1])/(mass_a);
    colvar.myder[i_c][i][2] =  (mass*dd_dra[2])/(mass_a);
  }
  for(i=colvar.list[i_c][0];i<colvar.list[i_c][0]+colvar.list[i_c][1];i++) {
    iat = colvar.cvatoms[i_c][i];
    mass=mtd_data->mass[iat];   
    colvar.myder[i_c][i][0] =  (mass*dd_drb[0])/(mass_b);
    colvar.myder[i_c][i][1] =  (mass*dd_drb[1])/(mass_b);
    colvar.myder[i_c][i][2] =  (mass*dd_drb[2])/(mass_b);
  }
  for(i=colvar.list[i_c][0]+colvar.list[i_c][1];i<colvar.list[i_c][0]+colvar.list[i_c][1]+colvar.list[i_c][2];i++) {
    iat = colvar.cvatoms[i_c][i];
    mass=mtd_data->mass[iat];   
    colvar.myder[i_c][i][0] =  (mass*dd_drh[0])/(mass_h);
    colvar.myder[i_c][i][1] =  (mass*dd_drh[1])/(mass_h);
    colvar.myder[i_c][i][2] =  (mass*dd_drh[2])/(mass_h);
  }
  for(i=colvar.list[i_c][0]+colvar.list[i_c][1]+colvar.list[i_c][2];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    mass=mtd_data->mass[iat];   
    colvar.myder[i_c][i][0] =  (mass*dd_drg[0])/(mass_g);
    colvar.myder[i_c][i][1] =  (mass*dd_drg[1])/(mass_g);
    colvar.myder[i_c][i][2] =  (mass*dd_drg[2])/(mass_g);
  }

  colvar.ss0[i_c] = mod_diff;

}

void PREFIX proj_on_axis_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i,j, iat;
  real mod_ba, mod_ba_3, mod_ha, mass_a, mass_b, mass_h, dot, mass;
  rvec r_ba,r_ha, sum_a, sum_b, sum_h;


  mass_a = mass_b = mass_h = 0.;
  sum_a[0] = sum_a[1] = sum_a[2] = 0.;
  sum_b[0] = sum_b[1] = sum_b[2] = 0.;
  sum_h[0] = sum_h[1] = sum_h[2] = 0.;

  j=0;
  // added optimized weights for axis fitting
  for(i=0;i<colvar.list[i_c][0];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][4]){mass=colvar.ow_weight[i_c*3][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    sum_a[0] += mass*mtd_data->pos[iat][0];
    sum_a[1] += mass*mtd_data->pos[iat][1];
    sum_a[2] += mass*mtd_data->pos[iat][2];
    if(colvar.intpar[i_c][4]){mass_a +=1;}
    else {mass_a += mass;}
  }
  j=0;
  for(i=colvar.list[i_c][0];i<colvar.list[i_c][0]+colvar.list[i_c][1];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][5]){mass=colvar.ow_weight[i_c*3+1][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    sum_b[0] += mass*mtd_data->pos[iat][0];
    sum_b[1] += mass*mtd_data->pos[iat][1];
    sum_b[2] += mass*mtd_data->pos[iat][2];
    if(colvar.intpar[i_c][5]){mass_b +=1;}
    else {mass_b += mass;}
  }
  j=0;
  for(i=colvar.list[i_c][0]+colvar.list[i_c][1];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][6]){mass=colvar.ow_weight[i_c*3+2][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    sum_h[0] += mass*mtd_data->pos[iat][0];
    sum_h[1] += mass*mtd_data->pos[iat][1];
    sum_h[2] += mass*mtd_data->pos[iat][2];
    if(colvar.intpar[i_c][6]){mass_h +=1;}
    else {mass_h += mass;}
  }


  sum_a[0] /= mass_a; sum_a[1] /= mass_a; sum_a[2] /= mass_a;
  sum_b[0] /= mass_b; sum_b[1] /= mass_b; sum_b[2] /= mass_b;
  sum_h[0] /= mass_h; sum_h[1] /= mass_h; sum_h[2] /= mass_h;


  if(colvar.cell_pbc[i_c]){
    minimal_image(sum_b, sum_a, &mod_ba, r_ba);
  } else {
    r_ba[0] = sum_b[0]-sum_a[0];
    r_ba[1] = sum_b[1]-sum_a[1];
    r_ba[2] = sum_b[2]-sum_a[2];
    mod_ba  = sqrt(r_ba[0]*r_ba[0]+r_ba[1]*r_ba[1]+r_ba[2]*r_ba[2]);
  };

  if(colvar.cell_pbc[i_c]){
    minimal_image(sum_h, sum_a, &mod_ha, r_ha);
  } else {
    r_ha[0] = sum_h[0]-sum_a[0];
    r_ha[1] = sum_h[1]-sum_a[1];
    r_ha[2] = sum_h[2]-sum_a[2];
  };

  // dot
  dot=r_ha[0]*r_ba[0]+r_ha[1]*r_ba[1]+r_ha[2]*r_ba[2];
  mod_ba_3=pow(mod_ba,3);

  // derivatives 
  j=0;
  for(i=0;i<colvar.list[i_c][0];i++) { 
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][4]){mass=colvar.ow_weight[i_c*3][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    colvar.myder[i_c][i][0] = mass/mass_a*(-(r_ba[0]+r_ha[0])/mod_ba+r_ba[0]*dot/mod_ba_3);  
    colvar.myder[i_c][i][1] = mass/mass_a*(-(r_ba[1]+r_ha[1])/mod_ba+r_ba[1]*dot/mod_ba_3); 
    colvar.myder[i_c][i][2] = mass/mass_a*(-(r_ba[2]+r_ha[2])/mod_ba+r_ba[2]*dot/mod_ba_3);
  }
  j=0;
  for(i=colvar.list[i_c][0];i<colvar.list[i_c][0]+colvar.list[i_c][1];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][5]){mass=colvar.ow_weight[i_c*3+1][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    colvar.myder[i_c][i][0] = mass/mass_b*(r_ha[0]/mod_ba-r_ba[0]*dot/mod_ba_3); 
    colvar.myder[i_c][i][1] = mass/mass_b*(r_ha[1]/mod_ba-r_ba[1]*dot/mod_ba_3); 
    colvar.myder[i_c][i][2] = mass/mass_b*(r_ha[2]/mod_ba-r_ba[2]*dot/mod_ba_3);
  }
  j=0;
  for(i=colvar.list[i_c][0]+colvar.list[i_c][1];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    if(colvar.intpar[i_c][6]){mass=colvar.ow_weight[i_c*3+2][j];j++;} 
    else{mass=mtd_data->mass[iat];}   
    colvar.myder[i_c][i][0] = mass/mass_h*r_ba[0]/mod_ba; 
    colvar.myder[i_c][i][1] = mass/mass_h*r_ba[1]/mod_ba; 
    colvar.myder[i_c][i][2] = mass/mass_h*r_ba[2]/mod_ba;
  }

  colvar.ss0[i_c] = dot/mod_ba; 

}
