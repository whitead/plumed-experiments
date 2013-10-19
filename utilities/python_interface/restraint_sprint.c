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

void PREFIX sprint_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int firstAtom, secondAtom, i, j, k, ix;
  int nn = colvar.nn[i_c], mm = colvar.mm[i_c], nat = colvar.natoms[i_c];
  int ind_atom = colvar.intpar[i_c][0];
  int n_distinct_masses = colvar.intpar[i_c][1];
  int list_atoms[nat], list_atoms2[nat], ind_sorted;
  rvec rij;
  real sqrtn=sqrt( (real) nat);
  real num, iden, mod_rij, rdist, func, dfunc, rNdist, rMdist;
  real ncoord, r_0[n_distinct_masses][n_distinct_masses], d_0[n_distinct_masses][n_distinct_masses];
  real threshold;
  threshold=pow(0.00001,1./(nn-mm));
  real cm_vec[nat*nat];
  real tmp1, tmp2, tmp3, tmpcv[nat], tmpmass[nat];
  int compute_cm;
  int element[nat];

// based on masses, label element: useful for pair-specific r_0 and d_0
  for(i=0;i<nat;i++){
    for(j=0;j<n_distinct_masses;j++){
      if(mtd_data->mass[colvar.cvatoms[i_c][i]] == colvar.realpar[i_c][j][0]) {
        element[i]=j; 
      }    
    }
  }
// fill matrix of pair-specific r_0 and d_0
  k=0;
  for(i=0;i<n_distinct_masses;i++){
    for(j=i;j<n_distinct_masses;j++){
      r_0[i][j]=colvar.realpar[i_c][k][1]; 
      r_0[j][i]=colvar.realpar[i_c][k][1]; 
      d_0[i][j]=colvar.realpar[i_c][k][2]; 
      d_0[j][i]=colvar.realpar[i_c][k][2]; 
      k++;
//      printf("||||| i j r_0 d_0 = %i %i %f %f\n",i,j,r_0[i][j],d_0[i][j]); // debug
    }
  }

  for (i=0;i<nat;i++) { for (ix=0;ix<3;ix++) { colvar.myder[i_c][i][ix]=0.; } }
    
  // check if cm and gradients have been already computed
  compute_cm=1;
  if(sprint_data.nat==nat){ // atoms number and index must match
    compute_cm=0;
    for (i=0;i<nat;i++) {
      if (colvar.cvatoms[i_c][i]!=colvar.cvatoms[sprint_data.icv][i]) { compute_cm=1; }
    }
    if (sprint_data.step!=colvar.it) { compute_cm=1; } // step must match
    if (logical.debug) { compute_cm=1; }
  }

                       ////////////////
  if (compute_cm==1) { // compute_cm //
                       ////////////////

//  printf("STEP %d CV %d : COMPUTING CM AND GRADIENTS!\n",colvar.it,i_c); // DEBUG
  // initialize backup variables
  if (sprint_data.lambda != NULL) { free_1dr_array_alloc(sprint_data.lambda); }
  if (sprint_data.cm     != NULL) { free_2dr_array_alloc(sprint_data.cm,sprint_data.nat); }
  if (sprint_data.grad   != NULL) { free_3dr_array_alloc(sprint_data.grad,sprint_data.nat,sprint_data.nat); }
  sprint_data.nat=nat;
  sprint_data.icv=i_c;
  sprint_data.step=colvar.it;
  sprint_data.lambda=float_1d_array_alloc(nat);
  sprint_data.cm    =float_2d_array_alloc(nat,nat);
  sprint_data.grad  =float_3d_array_alloc(nat,nat,3);

  // compute the contact matrix and gradients
  for (i=0;i<nat;i++) { sprint_data.cm[i][i]=0.; for (ix=0;ix<3;ix++) { sprint_data.grad[i][i][ix]=0.; } }
  for (i=0;i<nat-1;i++) { // loop over cm
    for (j=i+1;j<nat;j++) { // loop over cm
      // compute distance i-j
      firstAtom =colvar.cvatoms[i_c][i];
      secondAtom=colvar.cvatoms[i_c][j];
      if(colvar.cell_pbc[i_c]){
        minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
      } else {
        rij[0] = mtd_data->pos[firstAtom][0]-mtd_data->pos[secondAtom][0];
        rij[1] = mtd_data->pos[firstAtom][1]-mtd_data->pos[secondAtom][1];
        rij[2] = mtd_data->pos[firstAtom][2]-mtd_data->pos[secondAtom][2];
        mod_rij = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      }
//      printf("___ i j ei ej r_0 d_0 = %4i %4i %4i %4i %8.3f %8.3f\n",i,j,element[i],element[j],r_0[element[i]][element[j]],d_0[element[i]][element[j]]);  // debug
      // compute switching function
      rdist = (mod_rij-d_0[element[i]][element[j]])/r_0[element[i]][element[j]];
      /* analitic limit of the switching function */
      if(rdist<=0.){
        ncoord=1.;
        dfunc=0.;
      }else if(rdist>0.999999 && rdist<1.000001){
        ncoord=nn/mm;
        dfunc=0.5*nn*(nn-mm)/mm;
      }else if(rdist>threshold){
        ncoord=0.;
        dfunc=0.;
      }else{
        rNdist = pow(rdist, nn-1);
        rMdist = pow(rdist, mm-1);
        num = 1.-rNdist*rdist;
        iden = 1./(1.-rMdist*rdist);
        func = num*iden;
        ncoord = func;
        dfunc = ((-nn*rNdist*iden)+(func*(iden*mm)*rMdist))/(mod_rij*r_0[element[i]][element[j]]);
      }
      sprint_data.cm[i][j]=ncoord; sprint_data.cm[j][i]=ncoord; // contact matrix 
      for(ix=0;ix<3;ix++) {
        sprint_data.grad[i][j][ix] =  dfunc*rij[ix]; // = d cm_ij / d R_i = -d cm_ij / d R_j
        sprint_data.grad[j][i][ix] = -dfunc*rij[ix]; // = d cm_ji / d R_j = -d cm_ji / d R_i
	// in restraint_coord would be:   colvar.myder[i_c][i][ix] += +dfunc*rij[ix];
	// in restraint_coord would be:   colvar.myder[i_c][j][ix] += -dfunc*rij[ix];
      }
    } // end of loop over cm
  } // end of loop over cm
 
  // diagonalize the cm using routine ql77 in file restraint_spath.c
  for(i=0;i<nat;i++){
    for(j=0;j<nat;j++){
      cm_vec[nat*i+j]=sprint_data.cm[i][j]; // put matrix in a vector ...
    }
  }
  ql77(nat,cm_vec,sprint_data.lambda); // eigenvalues are sorted in ascending order
  for(j=0;j<nat;j++){
    for(i=0;i<nat;i++){
      sprint_data.cm[i][j]=cm_vec[nat*j+i]; // ... and back to matrix: columns are eigenvectors
    }
  }

  // make positive the principal eigenvector
  for(i=0;i<nat;i++){ sprint_data.cm[i][nat-1]=fabs(sprint_data.cm[i][nat-1]); }

    ////////////////////
  } // end compute_cm //
    ////////////////////

//  if (compute_cm==0){printf("STEP %d CV %d : NOT COMPUTING CM AND GRADIENTS!\n",colvar.it,i_c);} // DEBUG
//  printf("xxx eigenvalues    : ");for(i=0;i<nat;i++)printf("%4.2f ",sprint_data.lambda[i]);printf("\n"); //DEBUG
//  printf("xxx max eigenvector: ");for(i=0;i<nat;i++)printf("%4.2f ",sprint_data.cm[i][nat-1]);printf("\n"); //DEBUG


  // the CV
  for(i=0;i<nat;i++){ tmpcv[i]=sprint_data.cm[i][nat-1]; list_atoms[i]=i; 
//debug printf("i v %4d %8.4f\n",i,tmpcv[i]); 
  }

  // sort all atoms based on CV, irrespective of element
  realquicksort(tmpcv,list_atoms,0,nat-1);
//  for(i=0;i<nat;i++){ printf("i listatoms tmpcv cv[listatoms] = %4d %4d %8.4f %8.4f\n",i,list_atoms[i],tmpcv[i],sprint_data.cm[list_atoms[i]][nat-1]); } //debug

  // sort based on mass (arrange in groups: one sorted group per element)
  if(n_distinct_masses>1){
    for (i=0;i<nat;i++){ tmpmass[i]=mtd_data->mass[colvar.cvatoms[i_c][list_atoms[i]]]; }
    k=-1;
    for (j=0;j<n_distinct_masses;j++){
      for(i=0;i<nat;i++){
        if(tmpmass[i]==colvar.realpar[i_c][j][0]){
          k++; list_atoms2[k]=list_atoms[i];      
//          printf("i listatoms2 mass[listatoms2] cv[listatoms] = %4d %4d %8.4f %8.4f\n",k,list_atoms2[k],mtd_data->mass[colvar.cvatoms[i_c][list_atoms2[k]]],sprint_data.cm[list_atoms2[k]][nat-1]); //debug
        }
      }
    }
    for (i=0;i<nat;i++){ list_atoms[i]=list_atoms2[i]; }
  }

  // now the list is properly sorted
  ind_sorted=list_atoms[ind_atom];
//debug  printf("indatom indsorted = %4d %4d\n",ind_atom,ind_sorted);
//old colvar.ss0[i_c] = sprint_data.lambda[nat-1]*sqrtn*sprint_data.cm[ind_atom][nat-1]; 
  colvar.ss0[i_c] = sprint_data.lambda[nat-1]*sqrtn*sprint_data.cm[ind_sorted][nat-1]; 

  // the gradients
  for (i=0;i<nat-1;i++) { // loop over cm
    for (j=i+1;j<nat;j++) { // loop over cm
       tmp1=2.*sprint_data.cm[i][nat-1]*sprint_data.cm[j][nat-1]; 	// d lambda_N / d a_ij
       tmp2=0.; 			 	                // d v_ii^N / d a_ij
       for(k=0;k<nat-1;k++){
         tmp2+=sprint_data.cm[ind_sorted][k]*(sprint_data.cm[i][k]*sprint_data.cm[j][nat-1]+sprint_data.cm[j][k]*sprint_data.cm[i][nat-1])/(sprint_data.lambda[nat-1]-sprint_data.lambda[k]);
       }
       tmp3=tmp1*sprint_data.cm[ind_sorted][nat-1]+sprint_data.lambda[nat-1]*tmp2; // derivative of product lambda_N*v_ii^N
       for(ix=0;ix<3;ix++){
         colvar.myder[i_c][i][ix]+=sqrtn*tmp3*sprint_data.grad[i][j][ix];
         colvar.myder[i_c][j][ix]-=sqrtn*tmp3*sprint_data.grad[i][j][ix];
       }
     } // loop over cm
   } // loop over cm
}

//--------------------------------------------------------------------------------------------------------------

int PREFIX read_sprint(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw, iat, j, k, help, n_distinct_masses, npairs_r0, npairs_d0;
  double r_0, d_0;
  double delta = 0.0;
  real threshold, value;

  help=0;
  d_0=0.;
  npairs_d0=0;

  colvar.cell_pbc[count]=1; // default is PBC

  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR SPRINT\n"); help=1;}

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }
  iw=seek_word(word,"NN");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.nn[count]); } else { fprintf(fplog,"|- NEEDED NN KEYWORD FOR SPRINT\n"); help=1;}
  iw=seek_word(word,"MM");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.mm[count]);} else { fprintf(fplog,"|- NEEDED MM KEYWORD FOR SPRINT\n"); help=1;}

  iw=seek_word(word,"R_0");
  if(iw>=0) { 
    sscanf(word[iw+1],"%i", &npairs_r0); // number of pairs of elements (e.g. for hydrocarbons =3: C-C, C-H, H-H)
    if(npairs_r0>10){ fprintf(fplog,"ERROR: TOO MANY ELEMENT PAIRS !!!\n"); help=1; }
    for(i=0;i<npairs_r0;i++){
      // order: for atoms A B C -> A-A, A-B, A-C, B-B, B-C, C-C
      sscanf(word[iw+2+i],"%lf", &r_0);
      colvar.realpar[count][i][1] = (real) r_0;
    }
  } else { fprintf(fplog,"|- NEEDED R_0 KEYWORD FOR SPRINT\n"); help=1;}

  iw=seek_word(word,"D_0");
  if(iw>=0) { 
    sscanf(word[iw+1],"%i", &npairs_d0); // number of pairs of elements (e.g. for hydrocarbons =3: C-C, C-H, H-H)
    if(npairs_d0!=npairs_r0){ fprintf(fplog,"ERROR: ELEMENT PAIRS MUST BE THE SAME FOR R_0 AND D_0 !!!\n"); help=1; }
    for(i=0;i<npairs_d0;i++){
      // order: for atoms A B C -> A-A, A-B, A-C, B-B, B-C, C-C
      sscanf(word[iw+2+i],"%lf", &d_0);
      colvar.realpar[count][i][2] = (real) d_0;
    }
  }

  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}
  iw=seek_word(word,"INDEX"); // which atom to bias
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.intpar[count][0]); } 
  else { fprintf(fplog,"|- ERROR: NEEDED INDEX KEYWORD FOR SPRINT\n"); help=1;}
  if (colvar.intpar[count][0]<1||colvar.intpar[count][0]>colvar.natoms[count]) { 
    fprintf(fplog,"|- ERROR: THE ATOM INDEX MUST BE BETWEEN 1 AND %i\n",colvar.natoms[count]); help=1; }

  if(help){
          fprintf(fplog, "\n- SPRINT CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "SPRINT LIST <g1> NN 8 MM 16 R_0 1 0.7 SIGMA 1.0 INDEX 1\n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "            1 2 3 4 5 6 7 8 9 10\n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "       \n");
          fprintf(stderr, "PluMed dead with errors: check log file  \n");
          EXIT(); 
  }

//  colvar.r_0[count]      = (real) r_0;
//  colvar.d_0[count]      = (real) d_0;
  colvar.type_s[count]   = 55;

  fprintf(fplog, "%1i-SPRINT; ATOMS INVOLVED: %i; ", count+1, colvar.natoms[count]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON");
  else                       fprintf(fplog, " PBC OFF");
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  fprintf(fplog, "|--FUNCTIONAL FORM: (1-((dist_mat_rmsd-d_0)/r_0)^n) / (1-((dist_mat_rmsd-d_0)/r_0)^m) \n");
  fprintf(fplog, "|--PARAMETERS: n= %i m= %i index=%i\n", colvar.nn[count], colvar.mm[count], colvar.intpar[count][0]);
  fprintf(fplog, "|--  R_0 = "); for(i=0;i<npairs_r0;i++){ fprintf(fplog, "%f ",colvar.realpar[count][i][1]); }; fprintf(fplog, "\n");
  fprintf(fplog, "|--  D_0 = "); for(i=0;i<npairs_r0;i++){ fprintf(fplog, "%f ",colvar.realpar[count][i][2]); }; fprintf(fplog, "\n");
  threshold=pow(0.00001,1./(colvar.nn[count]-colvar.mm[count]));
  value=(1.-pow(threshold,colvar.nn[count]))/(1.-pow(threshold,colvar.mm[count]));
  fprintf(fplog, "|--CUTOFF VALUE: %f\n",value);
  // find largest r_0, d_0
  r_0=0.; d_0=0.;
  for(i=0;i<npairs_r0;i++){
    if(colvar.realpar[count][i][1]>r_0){ r_0=colvar.realpar[count][i][1]; }
    if(colvar.realpar[count][i][2]>d_0){ d_0=colvar.realpar[count][i][2]; }
  }
  fprintf(fplog, "|--CUTOFF DISTANCE: %f (LARGEST R_0 D_0 = %f %f)\n",threshold*r_0+d_0,r_0,d_0);

  colvar.intpar[count][0]--; // in C arrays start from zero...

  iat=0;
  fprintf(fplog,"|- SET MEMBERS: ");
  for(i=0;i<colvar.natoms[count];i++){fprintf(fplog," %d ",colvar.cvatoms[count][i]+1);if((i+1)%20==0)fprintf(fplog,"\n               ");}fprintf(fplog,"\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  sprint_data.nat=0; // useful later to know wether to compute cm and gradients...

  fprintf(fplog,"|- ATOM MASSES: ");
  for(i=0;i<colvar.natoms[count];i++){
    fprintf(fplog, "%5d %9.4f  ",colvar.cvatoms[count][i]+1,mtd_data.mass[colvar.cvatoms[count][i]]);
  } fprintf(fplog,"\n");
  // store distinct atomic masses in the vector colvar.realpar[count][:][0]
  n_distinct_masses=1;
  colvar.realpar[count][0][0]=mtd_data.mass[colvar.cvatoms[count][0]];
  for(i=1;i<colvar.natoms[count];i++){
    k=0;
    for(j=0;j<n_distinct_masses;j++){
      if(mtd_data.mass[colvar.cvatoms[count][i]] != colvar.realpar[count][j][0]){
        k++;
      }
    }
    if(k==n_distinct_masses) {
      n_distinct_masses++; 
      colvar.realpar[count][n_distinct_masses-1][0]=mtd_data.mass[colvar.cvatoms[count][i]];
    }
  }
  colvar.intpar[count][1]=n_distinct_masses;
// check that enough r_0, d_0 are given
  if(npairs_r0 != (n_distinct_masses*(n_distinct_masses+1))/2) {
    fprintf(fplog, "ERROR: IF N IS THE NUMBER OF DISTINCT ELEMENTS (MASSES) THEN YOU NEED N*(N+1)/2 R_0 PAIRS!!!\n");
    EXIT(); 
  }
  fprintf(fplog,"|- DISTINCT ATOM MASSES (NUMBER, VALUES): %2d ",n_distinct_masses);
  for(i=0;i<n_distinct_masses;i++){
    fprintf(fplog, "%9.4f  ",colvar.realpar[count][i][0]);
  } fprintf(fplog,"\n\n");

  return colvar.natoms[count]; 
}
