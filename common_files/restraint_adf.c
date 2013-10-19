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

#define ADF_CUTOFF 0.0001

void PREFIX adf_restraint(int i_c, int ignore_repeats, struct mtd_data_s *mtd_data) {
  int i, j, k, n, m, iat, jat, kat, nn, mm, n_c, ix;
  real lowB, upperB, dertmp, rdist, rNdist, rMdist, num, iden;
  real dot, theta, dtheta, weight;
  real mod_rij, mod_rij2, f_ij, df_ij, mod_rik, mod_rik2, f_ik, df_ik, r_0, d_0; 
  rvec rij, rik; real dijk[3][3], dij[2][3], dik[2][3];
  int cutflag[ colvar.list[i_c][2] ];

  // We calculate all beads for a given adf in the first call to this routine
  // so any further calls we just return without doing anything 
  for(i=0;i<i_c;i++){
     if( colvar.type_s[i]==52 && colvar.rdflab[i]==colvar.rdflab[i_c] ){
       if( ignore_repeats==1 ){ return; }   // This flag ensures we calculate each cv separately when debugging derivatives
     }
  }

  // Get all the parameters and store them in local variables
  nn = colvar.nn[i_c]; mm = colvar.mm[i_c];
  r_0 = colvar.r_0[i_c]; d_0 = colvar.d_0[i_c];

  // Initialize all colvars and derivatives to zero
  for(k=0;k<colvar.nconst;k++){
     if( colvar.type_s[k]==52 && colvar.rdflab[k]==colvar.rdflab[i_c] ){
        colvar.ss0[k]=0.0;
        for(i=0;i<colvar.natoms[k];i++){ 
          colvar.myder[k][i][0]=colvar.myder[k][i][1]=colvar.myder[k][i][2]=0.0; 
        }
     }  
  }  

  // Accumulate n(\theta) and its derivative
  for(i=0;i<colvar.list[i_c][0];i++){  // Loop over central atoms
     iat = colvar.cvatoms[i_c][i];

     // Loop over second atoms to determine which pairs of atoms can be cutoff
     n=0;
     for(k=colvar.list[i_c][0]+colvar.list[i_c][1];k<colvar.natoms[i_c];k++){
        kat = colvar.cvatoms[i_c][k];

        // Compute the distance between atom i and atom k and the vector connecting the two 
        if(colvar.cell_pbc[i_c]){ minimal_image(mtd_data->pos[iat], mtd_data->pos[kat], &mod_rik, rik); }
        else{
          rik[0] = mtd_data->pos[iat][0]-mtd_data->pos[kat][0];
          rik[1] = mtd_data->pos[iat][1]-mtd_data->pos[kat][1];
          rik[2] = mtd_data->pos[iat][2]-mtd_data->pos[kat][2];
          mod_rik  = sqrt(rik[0]*rik[0]+rik[1]*rik[1]+rik[2]*rik[2]);
        }

        // Is value of switching function very small
        cutflag[n]=0;
        rdist = (mod_rik-d_0)/r_0;
        if ( rdist>0. ) {  
           f_ik = ( 1 - pow(rdist,nn) ) / ( 1 - pow(rdist,mm) );
           if( f_ik<ADF_CUTOFF ){ cutflag[n]=1; }
        }
        n++;
     }

     // Now actually compute all the cvs
     m=0;
     for(j=colvar.list[i_c][0];j<colvar.list[i_c][0]+colvar.list[i_c][1];j++){   // Loop over first distance
        jat = colvar.cvatoms[i_c][j];

        // Cycle if two atoms are the same
        if( iat==jat ){ continue; }

        // Cycle if we don't need to include this vector (N.B. rdfNorm for adf's indicates whether or not lists 2 and 3 are identical)
        if( colvar.rdfNorm[i_c]==1 && cutflag[m]==1 ){ continue; } m++;     

        // Compute the distance between atom i and atom j and the vector connecting the two 
        if(colvar.cell_pbc[i_c]){ minimal_image(mtd_data->pos[iat], mtd_data->pos[jat], &mod_rij, rij); }
        else{
          rij[0] = mtd_data->pos[iat][0]-mtd_data->pos[jat][0];
          rij[1] = mtd_data->pos[iat][1]-mtd_data->pos[jat][1];
          rij[2] = mtd_data->pos[iat][2]-mtd_data->pos[jat][2];
          mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
        }
        mod_rij2 = mod_rij*mod_rij;    // Square of this distance

        // Compute the value of the switching function
        rdist = (mod_rij-d_0)/r_0;
        if(rdist<=0.){
           f_ij = 1; df_ij = 0.; 
        } else if (rdist>0.999999 && rdist<1.000001){
           f_ij = nn/mm; df_ij = 0.5*nn*(nn-mm)/mm;
        } else {
           rNdist = pow(rdist, nn-1); rMdist = pow(rdist, mm-1);
           num = 1.-rNdist*rdist;
           iden = 1./(1.-rMdist*rdist);
           f_ij = num*iden;
           df_ij = ((-nn*rNdist*iden)+(f_ij*(iden*mm)*rMdist))/(mod_rij*colvar.r_0[i_c]);
        }

        if( f_ij<ADF_CUTOFF ) { continue; }

        // Now compute the derivative of the switching function wrt the positions
        for(ix=0;ix<3;ix++) { dij[0][ix] = df_ij*rij[ix]; dij[1][ix] = -df_ij*rij[ix]; }

        n=0;
        for(k=colvar.list[i_c][0]+colvar.list[i_c][1];k<colvar.natoms[i_c];k++){ // Loop over the second distance
           kat = colvar.cvatoms[i_c][k];

           // Cutflag here ensures that we don't include very far away angles in the distribution
           // See code above
           if( cutflag[n]==1 || iat==kat || jat==kat ){ continue; } n++;

           // Compute the distance between atom i and atom k and the vector connecting the two 
           if(colvar.cell_pbc[i_c]){ minimal_image(mtd_data->pos[iat], mtd_data->pos[kat], &mod_rik, rik); }
           else{
             rik[0] = mtd_data->pos[iat][0]-mtd_data->pos[kat][0];
             rik[1] = mtd_data->pos[iat][1]-mtd_data->pos[kat][1];
             rik[2] = mtd_data->pos[iat][2]-mtd_data->pos[kat][2];
             mod_rik  = sqrt(rik[0]*rik[0]+rik[1]*rik[1]+rik[2]*rik[2]);
           }
           mod_rik2 = mod_rik*mod_rik;  // Square of this distance

           // Compute the value of the switching function
           rdist = (mod_rik-d_0)/r_0;
           if(rdist<=0.){
              f_ij = 1; df_ij = 0.;
           } else if (rdist>0.999999 && rdist<1.000001){
              f_ij = nn/mm; df_ij = 0.5*nn*(nn-mm)/mm;
           } else {
              rNdist = pow(rdist, nn-1); rMdist = pow(rdist, mm-1);
              num = 1.-rNdist*rdist;
              iden = 1./(1.-rMdist*rdist);
              f_ik = num*iden;
              df_ik = ((-nn*rNdist*iden)+(f_ik*(iden*mm)*rMdist))/(mod_rik*colvar.r_0[i_c]);
           }

           // Now compute the derivative of the switching function wrt the positions
           for(ix=0;ix<3;ix++) { dik[0][ix] = df_ik*rik[ix]; dik[1][ix] = -df_ik*rik[ix]; }
           
           // Compute this angle
           dot = iprod(rij,rik);
           num = 1 / ( mod_rij*mod_rik );
           iden = dot*num;     // The cosine of the angle
           theta = acos( iden );  // The angle
           dtheta = -num / sqrt( 1. - iden*iden );  // Derivative of the angle

           // The derivative of the angle with respect to the posisitions
           for(ix=0;ix<3;ix++) {
              dijk[0][ix] =  dtheta*( -dot/mod_rij2*rij[ix] - dot/mod_rik2*rik[ix] + rij[ix]+rik[ix] );
              dijk[1][ix] = -dtheta*( -dot/mod_rij2*rij[ix] + rik[ix] );
              dijk[2][ix] =  dtheta*(  dot/mod_rik2*rik[ix] - rij[ix] );
           }

           // Compute the weight of this term
           weight = f_ij*f_ik;

           // Now compute all the cvs
           for(n_c=0;n_c<colvar.nconst;n_c++){
              // Cycle if not a rdf bead of this type
              if( colvar.type_s[n_c]!=52 || colvar.rdflab[n_c]!=colvar.rdflab[i_c] ){ continue; }

              lowB = ( colvar.rdfBeadLower[n_c] - theta ) / ( sqrt(2.0) * colvar.rdfBeadWidth[n_c] ) ;
              upperB = ( colvar.rdfBeadUpper[n_c] - theta ) / ( sqrt(2.0) * colvar.rdfBeadWidth[n_c] ) ;  
              dertmp = ( exp( -lowB*lowB ) - exp( -upperB*upperB ) ) / ( sqrt(2.0*M_PI)*colvar.rdfBeadWidth[n_c] );

              iden = 0.5*( erf( upperB ) - erf( lowB ) );
              colvar.ss0[n_c] += weight*iden; 

              // The derivative of the distribution with respect to the atoms
              for(ix=0;ix<3;ix++) {
                 colvar.myder[n_c][i][ix] += weight*dertmp*dijk[0][ix] + iden*f_ik*dij[0][ix] + iden*f_ij*dik[0][ix];
                 colvar.myder[n_c][j][ix] += weight*dertmp*dijk[1][ix] + iden*f_ik*dij[1][ix];
                 colvar.myder[n_c][k][ix] += weight*dertmp*dijk[2][ix] + iden*f_ij*dik[1][ix];
              }
           }
        }
     }
  }

  // Normalize by dividing by the number of central atoms
  for(n_c=0;n_c<colvar.nconst;n_c++){
     if( colvar.type_s[n_c]!=52 || colvar.rdflab[n_c]!=colvar.rdflab[i_c] ){ continue; } 
     colvar.ss0[n_c] /= (double) colvar.list[i_c][0];
     for(i=0;i<colvar.natoms[n_c];i++){
         colvar.myder[n_c][i][0] /= (double) colvar.list[i_c][0];
         colvar.myder[n_c][i][1] /= (double) colvar.list[i_c][0];
         colvar.myder[n_c][i][2] /= (double) colvar.list[i_c][0]; 
     } 
  }

}

int PREFIX read_adf( char **word, int count, t_plumed_input *input, FILE *fplog ) {
  int i, iw, j, iat, k, help, rdflab, checkLists, natoms1, natoms2; 
  double delta = 0.0, upBound, lowBound, width, r_0, d_0;

  help=0; d_0=0;

  // Deal with periodic boundary conditions
  colvar.cell_pbc[count]=1; // default is PBC

  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;} 

  // This is the label for which rdf this bead is from
  iw = seek_word(word,"RDF_LABEL");
  if(iw>=0){ sscanf(word[iw+1],"%d",&rdflab); colvar.rdflab[count]=rdflab; }
  else{ fprintf(fplog,"|- NEEDED RDF_LABEL KEYWORD FOR ANGLE RDF\n"); help=1; }

  // Now we check whether there are other beads from this rdf that have been read in
  checkLists=-1; for(i=0;i<count;i++){ if( colvar.type_s[i]==52 && colvar.rdflab[i]==rdflab ){ checkLists=i; break; } }

  // Get atoms involved in the rdf
  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;
             j=plumed_get_group(word[iw+3],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][2]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR ADF\n"); help=1; }

  // Get the parameters for the switching function
  iw=seek_word(word,"NN");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.nn[count]); } else { fprintf(fplog,"|- NEEDED NN KEYWORD FOR ADF\n"); help=1;}
  iw=seek_word(word,"MM");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.mm[count]);} else { fprintf(fplog,"|- NEEDED MM KEYWORD FOR ADF\n"); help=1;}
  iw=seek_word(word,"R_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &r_0); } else { fprintf(fplog,"|- NEEDED R_0 KEYWORD FOR ADF\n"); help=1;}
  iw=seek_word(word,"D_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &d_0); }
   
  // Store the values
  colvar.r_0[count]      = (real) r_0; 
  colvar.d_0[count]      = (real) d_0; 

  // Check that the everything is consistent with other rdfs of this type
  if( checkLists!=-1){
    
    // Check the lists
    if( colvar.natoms[count]!=colvar.natoms[checkLists] ){
      fprintf(fplog,"|- MISLABELED ADF: ATOM LIST SIZES DO NOT MATCH\n"); help=1;  
    } else{
      for(i=0;i<colvar.natoms[count];i++){
         if( colvar.cvatoms[count][i]!=colvar.cvatoms[checkLists][i] ){
            fprintf(fplog,"|- MISLABELED ADF: ATOM LISTS DO NOT MATCH\n"); help=1;
         }
      }
    }
    
    // Now check the parameters for the switching functions 
    if( colvar.nn[count]!=colvar.nn[checkLists] ){
      fprintf(fplog,"|- MISLABELED ADF: NN VALUES DO NOT MATCH\n"); help=1; 
    }
    if( colvar.mm[count]!=colvar.mm[checkLists] ){
      fprintf(fplog,"|- MISLABELED ADF: MM  VALUES DO NOT MATCH\n"); help=1; 
    }
    if( colvar.r_0[count]!=colvar.r_0[checkLists] ){
      fprintf(fplog,"|- MISLABELED ADF: R_0 VALUES DO NOT MATCH\n"); help=1; 
    }
    if( colvar.d_0[count]!=colvar.d_0[checkLists] ){
      fprintf(fplog,"|- MISLABELED ADF: D_0 VALUES DO NOT MATCH\n"); help=1; 
    }
  }

  // Check whether or not lists 2 and 3 are the same
  colvar.rdfNorm[count]=1;
  natoms2=colvar.list[count][1] - colvar.list[count][0];
  natoms1=colvar.list[count][1] - colvar.list[count][0];   
  if( natoms1==natoms2 ){
    for(i=0;i<natoms1;i++){ 
       j=colvar.cvatoms[count][colvar.list[count][0]+i]; k=colvar.cvatoms[count][colvar.list[count][1]+i];
       if( j!=k ){ colvar.rdfNorm[count]=0; break; }
    }
  } else{ colvar.rdfNorm[count]=0; }

  // Get the start and end of this particular bead
  iw = seek_word(word,"RANGE");
  if(iw>=0){ sscanf(word[iw+1],"%lf",&lowBound); sscanf(word[iw+2],"%lf", &upBound); 
     colvar.rdfBeadLower[count]=lowBound; colvar.rdfBeadUpper[count]=upBound;
  } else{ fprintf(fplog,"|- NEEDED RANGE KEYWORD FOR RDF\n"); help=1;}

  // Get the width of this bead
  iw = seek_word(word,"WIDTH");
  if(iw>=0){ sscanf(word[iw+1],"%lf",&width); 
     colvar.rdfBeadWidth[count]=width;
  } else{ fprintf(fplog,"|- NEEDED WIDTH KEYWORD FOR RDF\n"); help=1;}

  if(help){
     fprintf(fplog,"\n-ADF CV: WRONG SYNTAX\n");
     fprintf(fplog,"e.g.:    \n");
     fprintf(fplog,"ADF RDF_LABEL 1 LIST <gr1> <gr2> <gr2> RANGE 2.0 3.0 WIDTH 1.5 R_0 1.0 D_0 0.0 NN 6 MM 12 \n");
     fprintf(fplog,"ADF RDF_LABEL 1 LIST <gr1> <gr2> <gr2> RANGE 2.0 3.0 WIDTH 1.5 R_0 1.0 D_0 0.0 NN 6 MM 12 \n");
     fprintf(fplog,"gr1-> \n");
     fprintf(fplog,"1 2 3 4  \n");
     fprintf(fplog,"gr1<- \n");
     fprintf(fplog,"gr2-> \n");
     fprintf(fplog,"1 2 3 4  \n");
     fprintf(fplog,"gr2<- \n");
     fprintf(fplog,"gr3-> \n");
     fprintf(fplog,"1 2 3 4  \n");
     fprintf(fplog,"gr3<- \n");
     fprintf(fplog,"         \n");
     plumed_error("PluMed dead with errors: check log file");
  }

  // Resize derivatives array
  snew(colvar.myder[count], colvar.natoms[count]);

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  fprintf(fplog, "%1i-BEAD FROM DISTRIBUTION OF ANGLES ABOUT (1st SET: %i ATOMS) ANGLES ARE BETWEEN (2nd SET: %i ATOMS) AND (3rd SET: %i ATOMS); \n", count+1, colvar.list[count][0], colvar.list[count][1], colvar.list[count][2]);
  if( colvar.rdfNorm[count]==1 ){ fprintf(fplog, "|--LISTS 2 AND 3 ARE IDENTICAL\n"); }
  fprintf(fplog, "|--LOWER BOUND %f, UPPER BOUND %f, WIDTH %f \n", colvar.rdfBeadLower[count], colvar.rdfBeadUpper[count], colvar.rdfBeadWidth[count]);
  fprintf(fplog, "|--PARAMETERS FOR SWITCHING FUNCTION : R_0 %f, D_0 %f, NN %d, MM %d \n", colvar.r_0[count], colvar.d_0[count], colvar.nn[count], colvar.mm[count] );  

  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 3rd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][2];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  return colvar.natoms[count];
}
