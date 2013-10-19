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

void PREFIX histogram_restraint(int i_c, struct mtd_data_s *mtd_data) {
  int i, j, k, j_c; real lowB, upperB, dertmp;


  k=0; colvar.ss0[i_c]=0.0;
  for(i=0;i<colvar.histo_ncv[i_c];i++){
     j_c = colvar.histo_cvlist[i_c][i];

     lowB = ( colvar.histo_low[i_c] - colvar.ss0[j_c] ) / ( sqrt(2.0) * colvar.histo_width[i_c] ) ;
     upperB = ( colvar.histo_high[i_c] - colvar.ss0[j_c] ) / ( sqrt(2.0) * colvar.histo_width[i_c] ) ;
     dertmp = ( exp( -lowB*lowB ) - exp( -upperB*upperB ) ) / ( sqrt(2.0*M_PI)*colvar.histo_width[i_c] );

     // Accumulate this colvar
     colvar.ss0[i_c]+=0.5*( erf( upperB ) - erf( lowB ) );

     // And derivatives 
     for(j=0;j<colvar.natoms[j_c];j++) {
        colvar.myder[i_c][k][0]=dertmp*colvar.myder[j_c][j][0];
        colvar.myder[i_c][k][1]=dertmp*colvar.myder[j_c][j][1];
        colvar.myder[i_c][k][2]=dertmp*colvar.myder[j_c][j][2];
        k++;
     }
  }

}

int PREFIX read_histogram( char **word, int count, t_plumed_input *input, FILE *fplog ) {
  int i, j, k, iw, ncolvar, ngrid, nspline;  int* cvlist;
  double smooth, delta = 0.0; int help=0;

  ncolvar=0; cvlist=NULL;

  // This reads in all the CVs used for histogram collective coordinates
  iw = seek_word(word,"CV_LIST"); 
  if(iw>=0){  ncolvar=plumed_get_group(word[iw+1],&cvlist,0,input,fplog); 
             if( ncolvar>count ){ plumed_error("NUMBER OF CVS USED TO CALCULATE HISTOGRAM"); }
  } 
  else{ fprintf(fplog,"|- NEEDED CV_LIST KEYWORD FOR HISTOGRAM\n"); help=1;}

  colvar.histo_ncv[count]=ncolvar; 
  snew( colvar.histo_cvlist[count], colvar.histo_ncv[count] );
  for(i=0;i<colvar.histo_ncv[count];++i){ colvar.histo_cvlist[count][i]=cvlist[i]; }

  // This is the sanity check for CV_LIST ( part 2 : check all the cvs used to construct histogram cvs are declared before any BESPOKE cvs)
  for(i=0;i<colvar.histo_ncv[count];i++){
     if( colvar.histo_cvlist[count][i]>count ){fprintf(fplog,"HISTOGRAM COMMAND MUST COME AFTER ALL CVS USED TO CREATE HISTOGRAM CV \n"); help=1; }
     // Number of atoms involved in this CV (sum of number of atoms involved in all other cvs)
     colvar.natoms[count]+=colvar.natoms[ colvar.histo_cvlist[count][i] ];
  }

  iw = seek_word(word,"RANGE");
  if(iw>=0) { 
     sscanf(word[iw+1], "%lf", &delta ); colvar.histo_low[count] = (real) delta; 
     sscanf(word[iw+2], "%lf", &delta ); colvar.histo_high[count] = (real) delta; 
  } 
  else{ fprintf(fplog,"|- NEEDED RANGE KEYWORD FOR HISTOGRAM\n"); help=1; }

  iw = seek_word(word,"WIDTH");
  if(iw>=0) { sscanf(word[iw+1], "%lf", &delta ); colvar.histo_width[count] = (real) delta; }
  else{ fprintf(fplog,"|- NEEDED WIDTH KEYWORD FOR HISTOGRAM\n"); help=1; }

  // Setup array of atoms involved in this cv
  snew( colvar.cvatoms[count], colvar.natoms[count] ); k=0;
  for(i=0;i<colvar.histo_ncv[count];i++){
     for(j=0;j<colvar.natoms[ colvar.histo_cvlist[count][i] ];j++){
        colvar.cvatoms[count][k]=colvar.cvatoms[ colvar.histo_cvlist[count][i] ][j];
        k++;
     }
  }
  // Setup derivatives array
  snew( colvar.myder[count], colvar.natoms[count] );

  if(help){ 
     fprintf(fplog,"\n-HISTOGRAM CV: WRONG SYNTAX\n");
     fprintf(fplog,"e.g.:    \n");
     fprintf(fplog,"TORSION LIST <g1> <g2> <g3> <g4>    \n");
     fprintf(fplog,"TORSION LIST <g2> <g3> <g4> <g5>    \n");
     fprintf(fplog,"TORSION LIST <g3> <g4> <g5> <g6>    \n");
     fprintf(fplog,"TORSION LIST <g4> <g5> <g6> <g7>    \n");
     fprintf(fplog,"HISTOGRAM CV_LIST <cv_list> RANGE 2.0 3.0 WIDTH 0.5\n");
     fprintf(fplog,"cvlist-> \n");
     fprintf(fplog,"1 2 3 4  \n");
     fprintf(fplog,"cvlist<- \n");
     fprintf(fplog,"         \n"); 
     plumed_error("PluMed dead with errors: check log file");
  }

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  fprintf(fplog, "%1i-HISTOGRAM COLLECTIVE COORDINATE CREATED FROM %i OTHER CVS \n", count+1, colvar.histo_ncv[count]);
  fprintf(fplog, "|--PARAMETERS: RANGE OF DESIRED VALUE = %f %f.  WIDTH OF GAUSSIAN FOR SMOOTHING %f. \n", colvar.histo_low[count], colvar.histo_high[count],colvar.histo_width[count] );
  fprintf(fplog,"|- COLVARS USED IN HISTOGRAM CV:");
  for(i=0;i<colvar.histo_ncv[count];i++){fprintf(fplog," %d ",colvar.histo_cvlist[count][i]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  return colvar.natoms[count];
}

// A separate test derivatives routine

void PREFIX histogram_testderivatives(int j_c, struct mtd_data_s *mtd_data)
{ 
  real teststep = 0.00001, invstep = 100000. ; 
  double testforce, analder, val1, val2;
  int i_c, i_bc, ix, i, j, k, iat, it;

  for(i=0;i<colvar.natoms[j_c];i++) {   // Loop over all atoms in the CV 
      iat = colvar.cvatoms[j_c][i];
      for(ix=0;ix<3;ix++) {              // Loop over all coordinates 
          mtd_data->pos[iat][ix] += teststep;
          for(i_bc=0;i_bc<colvar.histo_ncv[j_c];i_bc++) { // Loop over all cvs used to construct histo cv
              i_c=colvar.histo_cvlist[j_c][i_bc];
              switch(colvar.type_s[i_c]){
                case 1: dist_restraint(i_c, mtd_data); break;
                case 2: mindist_restraint(i_c, mtd_data); break;
                case 3: coord_restraint(i_c, mtd_data); break;
                case 4: angle_restraint(i_c, mtd_data); break;
                case 5: torsion_restraint(i_c, mtd_data); break;
                case 6: alfabeta_restraint(i_c, mtd_data); break;
                case 7: hbonds_restraint(i_c, mtd_data); break;
                case 8: dipole_restraint(i_c, mtd_data); break;
                case 11: radgyr_restraint(i_c, mtd_data); break;
                case 16: dihcor_restraint(i_c, mtd_data); break;
                case 20: waterbridge_restraint(i_c, mtd_data); break;
                case 30: spath_restraint(i_c, mtd_data); break;
                case 31: zpath_restraint(i_c, mtd_data); break;
                case 32: position_restraint(i_c, mtd_data); break;
                case 33: elstpot_restraint(i_c, mtd_data); break;
                case 34: puckering_restraint(i_c, mtd_data); break;
                case 36: helix_restraint(i_c, mtd_data); break;
                case 37: alpharmsd_restraint(i_c, mtd_data); break;
                case 38: antibetarmsd_restraint(i_c, mtd_data); break;
                case 39: parabetarmsd_restraint(i_c, mtd_data); break;
                //case 40: camshift_restraint(i_c, mtd_data); break;
                case 45: cmap_restraint(i_c, mtd_data); break;
                case 47: rdf_restraint(i_c, 1, mtd_data); break;
                case 52: adf_restraint(i_c, 0, mtd_data); break;
              }
          }
          histogram_restraint( j_c, mtd_data ); k=0; val1=colvar.ss0[j_c];

          mtd_data->pos[iat][ix] += -2.*teststep;
          for(i_bc=0;i_bc<colvar.histo_ncv[j_c];i_bc++) {
              i_c=colvar.histo_cvlist[j_c][i_bc];
              switch(colvar.type_s[i_c]){
               case 1: dist_restraint(i_c, mtd_data); break;
               case 2: mindist_restraint(i_c, mtd_data); break;
               case 3: coord_restraint(i_c, mtd_data); break;
               case 4: angle_restraint(i_c, mtd_data); break;
               case 5: torsion_restraint(i_c, mtd_data); break;
               case 6: alfabeta_restraint(i_c, mtd_data); break;
               case 7: hbonds_restraint(i_c, mtd_data); break;
               case 8: dipole_restraint(i_c, mtd_data); break;
               case 11: radgyr_restraint(i_c, mtd_data); break;
               case 16: dihcor_restraint(i_c, mtd_data); break;
               case 20: waterbridge_restraint(i_c, mtd_data); break;
               case 30: spath_restraint(i_c, mtd_data); break;
               case 31: zpath_restraint(i_c, mtd_data); break;
               case 32: position_restraint(i_c, mtd_data); break;
               case 33: elstpot_restraint(i_c, mtd_data); break;
               case 34: puckering_restraint(i_c, mtd_data); break;
               case 36: helix_restraint(i_c, mtd_data); break;
               case 37: alpharmsd_restraint(i_c, mtd_data); break;
               case 38: antibetarmsd_restraint(i_c, mtd_data); break;
               case 39: parabetarmsd_restraint(i_c, mtd_data); break;
               //case 40: camshift_restraint(i_c, mtd_data); break;
               case 45: cmap_restraint(i_c, mtd_data); break;
               case 47: rdf_restraint(i_c, 1, mtd_data); break;
               case 52: adf_restraint(i_c, 0, mtd_data); break;
              }
          }
          histogram_restraint( j_c, mtd_data ); val2=colvar.ss0[j_c];

          testforce = 0.5*((val1*invstep)-(val2*invstep)); 

          mtd_data->pos[iat][ix] += teststep;
          for(i_bc=0;i_bc<colvar.histo_ncv[j_c];i_bc++) {
              i_c=colvar.histo_cvlist[j_c][i_bc];
              switch(colvar.type_s[i_c]){
                case 1: dist_restraint(i_c, mtd_data); break;
                case 2: mindist_restraint(i_c, mtd_data); break;
                case 3: coord_restraint(i_c, mtd_data); break;
                case 4: angle_restraint(i_c, mtd_data); break;
                case 5: torsion_restraint(i_c, mtd_data); break;
                case 6: alfabeta_restraint(i_c, mtd_data); break;
                case 7: hbonds_restraint(i_c, mtd_data); break;
                case 8: dipole_restraint(i_c, mtd_data); break;
                case 11: radgyr_restraint(i_c, mtd_data); break;
                case 16: dihcor_restraint(i_c, mtd_data); break;
                case 20: waterbridge_restraint(i_c, mtd_data); break;
                case 30: spath_restraint(i_c, mtd_data); break;
                case 31: zpath_restraint(i_c, mtd_data); break;
                case 32: position_restraint(i_c, mtd_data); break;
                case 33: elstpot_restraint(i_c, mtd_data); break;
                case 34: puckering_restraint(i_c, mtd_data); break;
                case 36: helix_restraint(i_c, mtd_data); break;
                case 37: alpharmsd_restraint(i_c, mtd_data); break;
                case 38: antibetarmsd_restraint(i_c, mtd_data); break;
                case 39: parabetarmsd_restraint(i_c, mtd_data); break;
                //case 40: camshift_restraint(i_c, mtd_data); break;
                case 45: cmap_restraint(i_c, mtd_data); break;
                case 47: rdf_restraint(i_c, 1, mtd_data); break;
                case 52: adf_restraint(i_c, 0, mtd_data); break;
              }
          }
          histogram_restraint( j_c, mtd_data );
          
          analder=0; 
          for(j=0;j<colvar.natoms[j_c];j++) {
             if(colvar.cvatoms[j_c][j]==iat) {
                analder+=colvar.myder[j_c][j][ix]; 
             }
          }
          printf("Force atom %5i[%i] ** analytic %15.10f ** numeric %15.10f *** DELTA %15.10f\n",
              iat, ix, analder, testforce, analder-testforce); 
      }
  }

}
