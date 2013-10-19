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

// These two routines are only used in here

// double triangleFunc( double x ){
//     double result; result=fabs(x);
//     if(result<1.0){return 1.0-result;}
//     return 0.0;
// }   
// 
// double integratedTriangleFunc( double a, double b ){
//     double ia, ib; 
// 
//     if (b<=-1. || a >=1.) return 0.;
// 
//     if( a>-1.0 ){ ia=a; }else{ ia=-1.0; } 
//     if( b<1.0 ){ ib=b; } else{ ib=1.0; }
// 
//     return (ib*(2.-fabs(ib))-ia*(2.-fabs(ia)))*0.5;
// }  

void PREFIX rdf_restraint(int i_c, int ignore_repeats, struct mtd_data_s *mtd_data) {
  int i, j, k, iat, jat; 
  real lowB, upperB, dertmp, mod_rij, norm_factor; rvec rij;

//  fprintf(mtd_data->fplog,"New call %d\n",i_c );
//  fprintf(mtd_data->fplog, "ATOM 1 %f %f %f \n", mtd_data->pos[colvar.cvatoms[i_c][0]][0], mtd_data->pos[colvar.cvatoms[i_c][0]][1], mtd_data->pos[colvar.cvatoms[i_c][0]][2] );
//  fprintf(mtd_data->fplog,"ATOM 2 %f %f %f \n", mtd_data->pos[colvar.cvatoms[i_c][1]][0], mtd_data->pos[colvar.cvatoms[i_c][1]][1], mtd_data->pos[colvar.cvatoms[i_c][1]][2] );

  // We calculate all beads for a given rdf in the first call to this routine
  // so any further calls we just return without doing anything
  for(i=0;i<i_c;i++){ 
     if( colvar.type_s[i]==47 && colvar.rdflab[i]==colvar.rdflab[i_c] ){ 
       if( ignore_repeats==1 ){ return; }   // This flag ensures we calculate each cv separately when debugging derivatives 
     } 
  }

  // Initialize all colvars and derivatives to zero
  for(k=0;k<colvar.nconst;k++){
     if( colvar.type_s[k]==47 && colvar.rdflab[k]==colvar.rdflab[i_c] ){ 
       colvar.ss0[k]=0.0; 
       for(i=0;i<colvar.natoms[k];i++){ 
           colvar.myder[k][i][0]=colvar.myder[k][i][1]=colvar.myder[k][i][2]=0.0; 
       }
     }  
  }

  // Accumulate n(r) and its derivative 
  for(i=0;i<colvar.list[i_c][0];i++){   // Loop over central atoms
     iat = colvar.cvatoms[i_c][i];
     for(j=colvar.list[i_c][0];j<colvar.natoms[i_c];j++){   // Loop over other atoms
         jat = colvar.cvatoms[i_c][j];
         
         // Cycle if the two atoms are the same
         if( iat==jat ){ continue; }
         
         // Compute the distance between atom i and atom j and the derivative
         if(colvar.cell_pbc[i_c]){ minimal_image(mtd_data->pos[iat], mtd_data->pos[jat], &mod_rij, rij); }
         else{
           rij[0] = mtd_data->pos[iat][0]-mtd_data->pos[jat][0];
           rij[1] = mtd_data->pos[iat][1]-mtd_data->pos[jat][1];
           rij[2] = mtd_data->pos[iat][2]-mtd_data->pos[jat][2];
           mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
         }

         // Now loop over all colvars of this rdf and add the rdf and derivative
         for(k=0;k<colvar.nconst;k++){
            // Cycle if not a rdf bead of this type
            if( colvar.type_s[k]!=47 || colvar.rdflab[k]!=colvar.rdflab[i_c] ){ continue; }

            //lowB=(colvar.rdfBeadLower[k]-mod_rij) / colvar.rdfBeadWidth[k];
            //upperB=(colvar.rdfBeadUpper[k]-mod_rij) / colvar.rdfBeadWidth[k];

            //dertmp = triangleFunc( lowB ) / colvar.rdfBeadWidth[k] - triangleFunc( upperB ) / colvar.rdfBeadWidth[k];
            //colvar.ss0[k]+=integratedTriangleFunc( lowB, upperB );

            lowB = ( colvar.rdfBeadLower[k] - mod_rij ) / ( sqrt(2.0) * colvar.rdfBeadWidth[k] ) ;
            upperB = ( colvar.rdfBeadUpper[k] - mod_rij ) / ( sqrt(2.0) * colvar.rdfBeadWidth[k] ) ;  
            dertmp = ( exp( -lowB*lowB ) - exp( -upperB*upperB ) ) / ( sqrt(2.0*M_PI)*colvar.rdfBeadWidth[k] );
            colvar.ss0[k] += 0.5*( erf( upperB ) - erf( lowB ) ); 

//            fprintf( mtd_data->fplog,"Hello %d %d %d : %f %f %f %f %f %f\n", k, iat, jat, mod_rij, lowB, upperB, triangleFunc( lowB ), triangleFunc( upperB ), integratedTriangleFunc( lowB, upperB ) );

            // Accumulate the derivative on each distance pair
            colvar.myder[k][i][0] += dertmp*( rij[0]/mod_rij ); 
            colvar.myder[k][i][1] += dertmp*( rij[1]/mod_rij );
            colvar.myder[k][i][2] += dertmp*( rij[2]/mod_rij ); 
            colvar.myder[k][j][0] -= dertmp*( rij[0]/mod_rij );
            colvar.myder[k][j][1] -= dertmp*( rij[1]/mod_rij );
            colvar.myder[k][j][2] -= dertmp*( rij[2]/mod_rij );
         }
     }
  }

  real constant; constant = ( 4.0 * M_PI * colvar.list[i_c][0] ) / 3.0;

  // Convert n(r) to g(r) (assumes unit density)
  for(k=0;k<colvar.nconst;k++){
     if( colvar.type_s[k]==47 && colvar.rdflab[k]==colvar.rdflab[i_c] && colvar.rdfNorm[i_c]==1 ){
        norm_factor = constant * ( pow(colvar.rdfBeadUpper[k],3) - pow(colvar.rdfBeadLower[k],3) );
        colvar.ss0[k] /= norm_factor;
        for(i=0;i<colvar.natoms[k];i++){
            colvar.myder[k][i][0]/=norm_factor;
            colvar.myder[k][i][1]/=norm_factor;
            colvar.myder[k][i][2]/=norm_factor;
        }
     }
     else if( colvar.type_s[k]==47 && colvar.rdflab[k]==colvar.rdflab[i_c] ){
        colvar.ss0[k] /= (double) colvar.list[i_c][0];
        for(i=0;i<colvar.natoms[k];i++){
            colvar.myder[k][i][0] /= (double) colvar.list[i_c][0];
            colvar.myder[k][i][1] /= (double) colvar.list[i_c][0];
            colvar.myder[k][i][2] /= (double) colvar.list[i_c][0];
        }
     }
  }

}

int PREFIX read_rdf( char **word, int count, t_plumed_input *input, FILE *fplog ) {
  int i, iw, j, iat, help, rdflab, checkLists; 
  double delta = 0.0, upBound, lowBound, width;

  help=0;

  // Deal with periodic boundary conditions
  colvar.cell_pbc[count]=1; // default is PBC

  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;} 

  // This is the label for which rdf this bead is from
  iw = seek_word(word,"RDF_LABEL");
  if(iw>=0){ sscanf(word[iw+1],"%d",&rdflab); colvar.rdflab[count]=rdflab; }
  else{ fprintf(fplog,"|- NEEDED RDF_LABEL KEYWORD FOR RDF\n"); help=1; }

  // Now we check whether there are other beads from this rdf that have been read in
  checkLists=-1; for(i=0;i<count;i++){ if( colvar.type_s[i]==47 && colvar.rdflab[i]==rdflab ){ checkLists=i; break; } }

  // Get atoms involved in the rdf
  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR RDF\n"); help=1; }

  // Check that the lists are consistent with other rdfs of this type
  if( checkLists!=-1){
    if( colvar.natoms[count]!=colvar.natoms[checkLists] ){
      fprintf(fplog,"|- MISLABELED RDF: ATOM LIST SIZES DO NOT MATCH\n"); help=1;  
    }
    else{
      for(i=0;i<colvar.natoms[count];i++){
         if( colvar.cvatoms[count][i]!=colvar.cvatoms[checkLists][i] ){
            fprintf(fplog,"|- MISLABELED RDF: ATOM LISTS DO NOT MATCH\n"); help=1;
         }
      }
    }
  }

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

  // Choose whether or not to include the 4/3 pi ( r(i+1)**3 0 r(i)**3 ) normalization constant
  // default is to not include this term.
  iw = seek_word(word,"NORMALIZE");
  if(iw>=0){ colvar.rdfNorm[count]=1; }
  else{ colvar.rdfNorm[count]=0; }

  if(help){
     fprintf(fplog,"\n-RDF CV: WRONG SYNTAX\n");
     fprintf(fplog,"e.g.:    \n");
     fprintf(fplog,"RDF RDF_LABEL 1 LIST <gr1> <gr2> RANGE 2.0 3.0 WIDTH 1.5 \n");
     fprintf(fplog,"RDF RDF_LABEL 1 LIST <gr1> <gr2> RANGE 2.0 3.0 WIDTH 1.5 \n");
     fprintf(fplog,"gr1-> \n");
     fprintf(fplog,"1 2 3 4  \n");
     fprintf(fplog,"gr1<- \n");
     fprintf(fplog,"gr2-> \n");
     fprintf(fplog,"1 2 3 4  \n");
     fprintf(fplog,"gr2<- \n");
     fprintf(fplog,"         \n");
     plumed_error("PluMed dead with errors: check log file");
  }

  // Resize derivatives array
  snew(colvar.myder[count], colvar.natoms[count]);

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  fprintf(fplog, "%1i-RDF BEAD OF (1st SET: %i ATOMS) WRT (2nd SET: %i ATOMS); \n", count+1, colvar.list[count][0], colvar.list[count][1]);
  fprintf(fplog, "|--LOWER BOUND %f, UPPER BOUND %f, WIDTH %f \n", colvar.rdfBeadLower[count], colvar.rdfBeadUpper[count], colvar.rdfBeadWidth[count]);
  
  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  return colvar.natoms[count];
}
