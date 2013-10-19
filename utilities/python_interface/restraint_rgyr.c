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

void PREFIX radgyr_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i, j, firstAtom;
  rvec rcom, rdiff, pos0;
  real totmass, Rg, mod_rij,tmp;
  double gyration_tensor[3][3]; //gyration tensor must be calculated with double precision 
  double trace;
  double transf[3][3];  //transformation matrix
  double principal_components[3];
  double prefactor[3];
  double det;
  int pc_index;
  real tX[3];
  
  totmass = Rg = trace = rcom[0] = rcom[1] = rcom[2] = 0.;

  /* Initialise gyration tensor matrix */
  for(i=0;i<3;i++)
   for(j=0;j<3;j++)
     gyration_tensor[i][j]=0.0;

  firstAtom = colvar.cvatoms[i_c][0];
  if (colvar.mm[i_c]==0) totmass += mtd_data->mass[firstAtom];   //mass-weighted
  else                   totmass += 1.0;                         //non-weighted
  for(i=0;i<3;i++) pos0[i] = mtd_data->pos[firstAtom][i];

  /*Calculate center of mass/ geometrical center*/
  for(i=1;i<colvar.natoms[i_c];i++){
    firstAtom = colvar.cvatoms[i_c][i];
    if (colvar.mm[i_c]==0) totmass += mtd_data->mass[firstAtom];
    else                   totmass +=1.0; 
    if(colvar.cell_pbc[i_c]) {
      minimal_image(mtd_data->pos[firstAtom], pos0, &mod_rij, rdiff);
      if (colvar.mm[i_c]==0) for(j=0;j<3;j++) rcom[j] += mtd_data->mass[firstAtom]*rdiff[j];
      else                   for(j=0;j<3;j++) rcom[j] += rdiff[j];
    } else {
      if (colvar.mm[i_c]==0)  for(j=0;j<3;j++) rcom[j] += mtd_data->mass[firstAtom]*(mtd_data->pos[firstAtom][j]-pos0[j]);
      else                    for(j=0;j<3;j++) rcom[j] += mtd_data->pos[firstAtom][j]-pos0[j];
    }
  }
  for(j=0;j<3;j++) rcom[j] = rcom[j]/totmass+pos0[j];
  
  /* Calculate Rg / gyration tensor */
  for(i=0;i<colvar.natoms[i_c];i++){
    firstAtom = colvar.cvatoms[i_c][i];
    if(colvar.cell_pbc[i_c]) {
     minimal_image(mtd_data->pos[firstAtom], rcom, &mod_rij, rdiff);
    } else {
     for(j=0;j<3;j++) rdiff[j] = mtd_data->pos[firstAtom][j]-rcom[j];
     mod_rij  = sqrt(rdiff[0]*rdiff[0]+rdiff[1]*rdiff[1]+rdiff[2]*rdiff[2]);    
    }
 
    if(colvar.nn[i_c]<2){    // calculate Rg directly, no tensor needed
      if (colvar.mm[i_c]==0){
             Rg += mtd_data->mass[firstAtom]*mod_rij*mod_rij;
             for(j=0;j<3;j++) colvar.myder[i_c][i][j] = rdiff[j]*mtd_data->mass[firstAtom]; 
      }
      else{
             Rg += mod_rij*mod_rij;
             for(j=0;j<3;j++) colvar.myder[i_c][i][j] = rdiff[j];       
      }
    }else{
      if (colvar.mm[i_c]==0){  //calculate gyration tensor
             gyration_tensor[0][0]+=mtd_data->mass[firstAtom]*rdiff[0]*rdiff[0];
             gyration_tensor[1][1]+=mtd_data->mass[firstAtom]*rdiff[1]*rdiff[1];
             gyration_tensor[2][2]+=mtd_data->mass[firstAtom]*rdiff[2]*rdiff[2];
             gyration_tensor[0][1]+=mtd_data->mass[firstAtom]*rdiff[0]*rdiff[1];
             gyration_tensor[0][2]+=mtd_data->mass[firstAtom]*rdiff[0]*rdiff[2];
             gyration_tensor[1][2]+=mtd_data->mass[firstAtom]*rdiff[1]*rdiff[2];
      }
      else{
             gyration_tensor[0][0]+=rdiff[0]*rdiff[0];
             gyration_tensor[1][1]+=rdiff[1]*rdiff[1];
             gyration_tensor[2][2]+=rdiff[2]*rdiff[2];
             gyration_tensor[0][1]+=rdiff[0]*rdiff[1];
             gyration_tensor[0][2]+=rdiff[0]*rdiff[2];
             gyration_tensor[1][2]+=rdiff[1]*rdiff[2];
      }
    }
  }
  
  /*Evaluate gyration- and inertia- tensor based CV*/
  if((colvar.nn[i_c]>1)&&(colvar.nn[i_c]<11)){
      //dsyevj3(gyration_tensor,transf,principal_components); //diagonalize gyration tensor
      rank3_to_ql77(gyration_tensor, transf, principal_components);
      //sort eigenvalues and eigenvectors
      if (principal_components[0]<principal_components[1]){
      tmp=principal_components[0]; principal_components[0]=principal_components[1]; principal_components[1]=tmp;
      for (i=0; i<3; i++){tmp=transf[i][0]; transf[i][0]=transf[i][1]; transf[i][1]=tmp;}
      }
      if (principal_components[1]<principal_components[2]){
      tmp=principal_components[1]; principal_components[1]=principal_components[2]; principal_components[2]=tmp;
      for (i=0; i<3; i++){tmp=transf[i][1]; transf[i][1]=transf[i][2]; transf[i][2]=tmp;}
      }
      if (principal_components[0]<principal_components[1]){
      tmp=principal_components[0]; principal_components[0]=principal_components[1]; principal_components[1]=tmp;
      for (i=0; i<3; i++){tmp=transf[i][0]; transf[i][0]=transf[i][1]; transf[i][1]=tmp;}      
      }
      //calculate determinant of transformation matrix
      det=transf[0][0]*transf[1][1]*transf[2][2]+transf[0][1]*transf[1][2]*transf[2][0]+transf[0][2]*transf[1][0]*transf[2][1]- transf[0][2]*transf[1][1]*transf[2][0]-transf[0][1]*transf[1][0]*transf[2][2]-transf[0][0]*transf[1][2]*transf[2][1];
      if (det<0) 
	  for (j=0;j<3;j++) transf[j][2]=-transf[j][2]; //trasformation matrix for rotation must have positive determinant, otherwise multiply one column by (-1)
      det=transf[0][0]*transf[1][1]*transf[2][2]+transf[0][1]*transf[1][2]*transf[2][0]+transf[0][2]*transf[1][0]*transf[2][1]- transf[0][2]*transf[1][1]*transf[2][0]-transf[0][1]*transf[1][0]*transf[2][2]-transf[0][0]*transf[1][2]*transf[2][1];	    
      if (fabs(det-1)>0.001) plumed_error("Plumed Error: Cannot diagonalize gyration tensor\n"); //check again, if det(transf)!=1 something is wrong, die
      prefactor[0]=prefactor[1]=prefactor[2]=0;
      
      if (colvar.nn[i_c]<5){        //GTPC_1, GTPC_2, GTPC_3 (S'_1, S'_2, S'_3 in paper )
            pc_index=colvar.nn[i_c]-2; //index of principal component
            colvar.ss0[i_c]=sqrt(principal_components[pc_index]/totmass);
            if (colvar.ss0[i_c]>1e-6) prefactor[pc_index]=1.0/(totmass*colvar.ss0[i_c]); //some parts of derivate
      }
      else switch(colvar.nn[i_c]){
	    case 8:        //the smallest principal radius of gyration
                   colvar.ss0[i_c]=sqrt((principal_components[1]+principal_components[2])/totmass);
		   if (colvar.ss0[i_c]>REAL_EPS){
                   prefactor[1]=1.0/(totmass*colvar.ss0[i_c]);
                   prefactor[2]=1.0/(totmass*colvar.ss0[i_c]);
		   }
                   break;
	    case 9:       //the midle principal radius of gyration
                   colvar.ss0[i_c]=sqrt((principal_components[0]+principal_components[2])/totmass);
		   if (colvar.ss0[i_c]>REAL_EPS){
                   prefactor[0]=1.0/(totmass*colvar.ss0[i_c]);
                   prefactor[2]=1.0/(totmass*colvar.ss0[i_c]);
		   }
                   break;
	    case 10:      //the largest principal radius of gyration
                   colvar.ss0[i_c]=sqrt((principal_components[0]+principal_components[1])/totmass);
		   if (colvar.ss0[i_c]>REAL_EPS){
                   prefactor[0]=1.0/(totmass*colvar.ss0[i_c]);
                   prefactor[1]=1.0/(totmass*colvar.ss0[i_c]);
		   }
                   break;             
	    case 5:       //ASPHERICITY (b')
                  colvar.ss0[i_c]=sqrt((principal_components[0]-0.5*(principal_components[1]+principal_components[2]))/totmass); 
	          if (colvar.ss0[i_c]*totmass>REAL_EPS){   //avoid division by zero 
                  prefactor[0]= 1.0/(totmass*colvar.ss0[i_c]);
                  prefactor[1]=-0.5/(totmass*colvar.ss0[i_c]);
                  prefactor[2]=-0.5/(totmass*colvar.ss0[i_c]);
		  }
		  break;
      
	    case 6:     //ACYLINDRICITY (c')
                  colvar.ss0[i_c]=sqrt((principal_components[1]-principal_components[2])/totmass); 
	          if (colvar.ss0[i_c]*totmass>REAL_EPS){   //avoid division by zero  
                  prefactor[1]= 1.0/(totmass*colvar.ss0[i_c]);
                  prefactor[2]=-1.0/(totmass*colvar.ss0[i_c]);
		  }
                  break;
	    case 7:    //KAPPA2 - relative shape anisotropy
                  trace=principal_components[0]+principal_components[1]+principal_components[2];
                  tmp=principal_components[0]*principal_components[1]+ principal_components[1]*principal_components[2]+ principal_components[0]*principal_components[2]; 
                  colvar.ss0[i_c]=1.0-3*(tmp/(trace*trace));
		  if (colvar.ss0[i_c]>REAL_EPS){
                  prefactor[0]= -3*((principal_components[1]+principal_components[2])-2*tmp/trace)/(trace*trace) *2;
                  prefactor[1]= -3*((principal_components[0]+principal_components[2])-2*tmp/trace)/(trace*trace) *2;
                  prefactor[2]= -3*((principal_components[0]+principal_components[1])-2*tmp/trace)/(trace*trace) *2;
		  }
                  break;
	    default: break;
      }      
      /* Calculate the final part of the gradient and project it back to the original coordinate frame */
      for(i=0;i<colvar.natoms[i_c];i++){
                   firstAtom = colvar.cvatoms[i_c][i];
                   if (colvar.cell_pbc[i_c]) {
                             minimal_image(mtd_data->pos[firstAtom], rcom, &mod_rij, rdiff);
                   } else {
                   for (j=0;j<3;j++) rdiff[j] = mtd_data->pos[firstAtom][j]-rcom[j];
                   }
                   for (j=0;j<3;j++) tX[j]=transf[0][j]*rdiff[0]+transf[1][j]*rdiff[1]+transf[2][j]*rdiff[2]; //project atomic postional vectors to diagonalized frame
		   /* Complete the gradient in diagonalized frame and back-project to original coordinate frame */  
                   if (colvar.mm[i_c]==0) for (j=0;j<3;j++) colvar.myder[i_c][i][j]=mtd_data->mass[firstAtom]*(prefactor[0]*transf[j][0]*tX[0]+prefactor[1]*transf[j][1]*tX[1]+prefactor[2]*transf[j][2]*tX[2]);
                   else                   for (j=0;j<3;j++) colvar.myder[i_c][i][j]=                          (prefactor[0]*transf[j][0]*tX[0]+prefactor[1]*transf[j][1]*tX[1]+prefactor[2]*transf[j][2]*tX[2]);
        }
  }

/* Trace of the gyration tensor */
    if(colvar.nn[i_c]==0){
      colvar.ss0[i_c] = 2*Rg;
      for(i=0;i<colvar.natoms[i_c];i++) {
       for(j=0;j<3;j++) colvar.myder[i_c][i][j] *= 4;  
      }
/* Gyration radius */
    }else if(colvar.nn[i_c]==1){
      colvar.ss0[i_c] = sqrt(Rg/totmass);
      for(i=0;i<colvar.natoms[i_c];i++) {
       for(j=0;j<3;j++) colvar.myder[i_c][i][j] /= colvar.ss0[i_c]*totmass; 
      }
    }
}

//------------------------------------------------------------------------------------------------

int PREFIX read_rgyr(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw , j;
  double delta = 0.0;
  char string[400];
  int help;

  help=0;
  colvar.nn[count] = 0;
  colvar.mm[count] = 0; //0 = mass-weighted, other means no-weighted

  iw=seek_word(word,"RGYR");
  if(iw>=0)colvar.nn[count] = 1; //default

  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR RGYR\n"); help=1;} 

  iw=seek_word(word,"TYPE");
  if(iw>0){
    if(seek_word2(word,"TRACE",iw)>iw) colvar.nn[count] = 0;
    if(seek_word2(word,"RGYR",iw)>iw)  colvar.nn[count] = 1;
    if(seek_word2(word,"GTPC_1",iw)>iw)  colvar.nn[count] = 2;
    if(seek_word2(word,"GTPC_2",iw)>iw)  colvar.nn[count] = 3;
    if(seek_word2(word,"GTPC_3",iw)>iw)  colvar.nn[count] = 4;
    if(seek_word2(word,"ASPHERICITY",iw)>iw) colvar.nn[count] = 5;    
    if(seek_word2(word,"ACYLINDRICITY",iw)>iw) colvar.nn[count] = 6;        
    if(seek_word2(word,"KAPPA2",iw)>iw) colvar.nn[count] = 7; 
    if(seek_word2(word,"RGYR_3",iw)>iw)  colvar.nn[count] = 8;               
    if(seek_word2(word,"RGYR_2",iw)>iw)  colvar.nn[count] = 9;
    if(seek_word2(word,"RGYR_1",iw)>iw)  colvar.nn[count] = 10;


  }

  iw=seek_word(word,"PBC"); 
  if(iw>0) colvar.cell_pbc[count] = 1;
 
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  iw=seek_word(word,"MASS-WEIGHTED"); 
  if(iw>0) colvar.mm[count] = 0;
  
  iw=seek_word(word,"NO-WEIGHTED"); 
  if(iw>0) colvar.mm[count] = 1;  

  if(help){
          printf("%i \n",help);
          fprintf(fplog, "\n-gyration/RGYR CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "gyration/RGYR LIST <g1> SIGMA 0.1 [PBC] [MASS-WEIGHTED (default) /NO-WEIGHTED]\n [TYPE RGYR (default)/TRACE/RGYR_1/RGYR_2/RGYR_3/GTPC_1/GTPC_2/GTPC_3/ASPHERICITY/ACYLINDRICITY/KAPPA2]\n");
          fprintf(fplog, "         g1->            \n");
          fprintf(fplog, "         6 10 16         \n");
          fprintf(fplog, "         LOOP 20 30 2    \n");
          fprintf(fplog, "         g1<-            \n");
          fprintf(fplog, "                         \n");
          plumed_error("PluMeD dead with errors: check log file");
  }

  colvar.type_s[count]   = 11;

  snew(colvar.myder[count], colvar.natoms[count]);

  if(colvar.nn[count] == 0)
    fprintf(fplog,"%i-TRACE OF THE GYRATION TENSOR; ", count+1); 

  if(colvar.nn[count] == 1)
    fprintf(fplog,"%i-GYRATION RADIUS (Rg); ", count+1); 
  
  if(colvar.nn[count] == 2)
    fprintf(fplog,"%i-THE LARGEST PRINCIPAL MOMENT OF THE GYRATION TENSOR (S'_1); ", count+1); 

  if(colvar.nn[count] == 3)
    fprintf(fplog,"%i-THE MIDDLE PRINCIPAL MOMENT OF THE GYRATION TENSOR (S'_2); ", count+1); 
  
  if(colvar.nn[count] == 4)
    fprintf(fplog,"%i-THE SMALLEST PRINCIPAL MOMENT OF THE GYRATION TENSOR (S'_3); ", count+1); 
    
  if(colvar.nn[count] == 5)
    fprintf(fplog,"%i-THE ASPHERICITY (b'); ", count+1); 

  if(colvar.nn[count] == 6)
    fprintf(fplog,"%i-THE ACYLINDRICITY (c'); ", count+1); 

  if(colvar.nn[count] == 7)
    fprintf(fplog,"%i-THE RELATIVE SHAPE ANISOTROPY (kappa^2); ", count+1); 

  if(colvar.nn[count] == 8)
    fprintf(fplog,"%i-THE SMALLEST PRINCIPAL RADIUS OF GYRATION (r_g3); ", count+1); 

  if(colvar.nn[count] == 9)
    fprintf(fplog,"%i-THE MIDDLE PRINCIPAL RADIUS OF GYRATION (r_g2); ", count+1);

  if(colvar.nn[count] == 10)
    fprintf(fplog,"%i-THE LARGEST PRINCIPAL RADIUS OF GYRATION (r_g1); ", count+1);

  fprintf(fplog,"ATOMS INVOLVED: %i; ",colvar.natoms[count]);

  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  
  if (colvar.cell_pbc[count]) {fprintf(fplog,"|- PBC is ON \n");}
  else {fprintf(fplog,"|- PBC is OFF \n");} 
  
  if (colvar.mm[count]==0) {fprintf(fplog,"|- MASS-WEIGHTED \n");}
  else {fprintf(fplog,"|- NO-WEIGHTED \n");}   

  fprintf(fplog,"|- SET MEMBERS: ");
  for(i=0;i<colvar.natoms[count];i++){fprintf(fplog," %d ",colvar.cvatoms[count][i]+1);if((i+1)%20==0)fprintf(fplog,"\n               ");}fprintf(fplog,"\n\n");

  return colvar.natoms[count];
}

void PREFIX rank3_to_ql77(double in[3][3], double evector[3][3], double evalue[3]){
 int i,j;
 double ll[9];
 double *evalue2;  
 evalue2=(double *)malloc(3*sizeof(double));

 for(i=0;i<3;i++) for(j=0;j<3;j++) ll[3*i+j]=in[i][j];

// call diagonalization
 ql77(3,ll,evalue2);

//back to square representation: columns have eigenvectors
 for(j=0;j<3;j++){
    evalue[j]= evalue2[j];
    for(i=0;i<3;i++) evector[i][j]= ll[3*j+i]; 
 }
 free(evalue2); 
}
