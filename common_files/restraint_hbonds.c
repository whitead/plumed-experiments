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

void PREFIX hbonds_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int j, ix, i, kk, firstAtom, secondAtom;
  double nn, mm;
  rvec rij;
  real mod_rij, rdist;
  double nhbond, r6dist, r12dist, num, iden;
  double rijfmod, deriv;
  int firstRes = 1, secondRes = -1;  // ### Default residue numbers always different.

  nhbond = 0.;
  nn = (double) colvar.nn[i_c];
  mm = (double) colvar.mm[i_c];
  for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] = 0.;

  for(i=0;i<colvar.list[i_c][0];i++) {						// cycle over acceptors
    firstAtom = colvar.cvatoms[i_c][i];						// acceptors atom
    if (colvar.list[i_c][2] == colvar.list[i_c][0]) {
    	firstRes = colvar.cvatoms[i_c][i + colvar.natoms[i_c]]; // #### acceptors atom
    }
    for(j=colvar.list[i_c][0];j<colvar.natoms[i_c];j++) {	                // cycle over hydrogens
      kk=i-j+colvar.list[i_c][0];
      if(colvar.type[i_c]==1) {
        if(abs(kk)<=4) continue;
        if(kk%2==0) continue;
      } else if(colvar.type[i_c]==2) {
        if(abs(kk)!=4) continue;
      } else if(colvar.type[i_c]==3) {
        if(abs(kk)<=4) continue;
        if(kk%2==1) continue;
      } else if(colvar.type[i_c]==4) {
        if(kk!=0) continue;   
      } else if(colvar.type[i_c]==6) {    // #### This type is undocumented, so we put it to six for consistency
        if(abs(kk)<=4) continue; 
      }
      secondAtom = colvar.cvatoms[i_c][j];			// hydrogen atom
      if (colvar.list[i_c][3] == colvar.list[i_c][1]) {
    	  secondRes = colvar.cvatoms[i_c][j + colvar.natoms[i_c] ];
      }

      if(firstAtom==secondAtom) continue;
      // ##### There must be no force for donor and acceptor atoms on the same peptide group.
      if(colvar.type[i_c]==5) {
          if (firstRes == secondRes) continue;
    	  if(abs(firstAtom - secondAtom)<=5) continue;
      }
      if(colvar.cell_pbc[i_c]){ // distance acceptor/donor
        minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
      } else {
        rij[0] = mtd_data->pos[firstAtom][0]-mtd_data->pos[secondAtom][0];
        rij[1] = mtd_data->pos[firstAtom][1]-mtd_data->pos[secondAtom][1];
        rij[2] = mtd_data->pos[firstAtom][2]-mtd_data->pos[secondAtom][2];
        mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      };
      rdist = mod_rij/colvar.r_0[i_c];				// weighted distance
      if(rdist>0.99998 && rdist<1.00002) rdist = 0.99999;	// to keep the function continuos
      r6dist = pow(rdist, nn);	
      num = 1.-r6dist;						// numerator
      r12dist = pow(rdist, mm);	
      iden = 1./(1.-r12dist);					// denominator

      nhbond += num*iden;
      for(ix=0;ix<3;ix++) {
        rijfmod = (double) rij[ix]/(mod_rij*mod_rij);
        deriv = -rijfmod*(-mm*r12dist+r6dist*(nn+(mm-nn)*r12dist))*iden*iden;
        colvar.myder[i_c][i][ix] += (real) deriv;
        colvar.myder[i_c][j][ix] += (real) -deriv;
      }
    }
  }
  colvar.ss0[i_c] = (real) nhbond;
  
}

// --------------------------------------------------------------------------------------------------------------

int PREFIX read_hbonds(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  double delta = 0.0;
  double r_0;
  int iat, i, j, iw, help;
  char string[400];

  help=0;
  colvar.nn[count] = 6;
  colvar.mm[count] = 12;

#if defined (PLUMED_GROMACS)
  r_0 = 0.25;
#else
  r_0 = 2.5;
#endif

  colvar.cell_pbc[count]=1; // default is PBC

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
  }  else {fprintf(fplog,"|- NEEDED LIST KEYWORD FOR HBONDS\n"); help=1;}
  // ##### List of residue numbers for excluding intra-residue h-bonds.
  iw = seek_word(word,"RESLIST");
  if(iw>=0){
               j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
               colvar.list[count][2]=j;
               if (j != colvar.list[count][0] ) plumed_error("HBONDS RESLIST : residue numbers should correspond to atoms!");
               j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count] + j ,input,fplog);
               colvar.list[count][3]=j;
               if (j != colvar.list[count][1] ) plumed_error("HBONDS RESLIST : residue numbers should correspond to atoms!");
  }else{
			  colvar.list[count][2]=0;
			  colvar.list[count][3]=0;
  }

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }
  iw=seek_word(word,"NN");
  if(iw>=0) {sscanf(word[iw+1],"%i", &colvar.nn[count]); }
  iw=seek_word(word,"MM");
  if(iw>=0) {sscanf(word[iw+1],"%i", &colvar.mm[count]); }
  iw=seek_word(word,"R_0");
  if(iw>=0) {sscanf(word[iw+1],"%lf", &r_0); } 
  iw=seek_word(word,"TYPE");
  if(iw>=0) {sscanf(word[iw+1],"%i", &colvar.type[count]); } else {fprintf(fplog,"|- NEEDED TYPE KEYWORD FOR HBONDS\n"); help=1;}

  if(help){
          fprintf(fplog, "\n-HBONDS CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "HBONDS LIST <H> <O> TYPE 0 SIGMA 0.1 [NN 6 MM 12 R_0 2.5]\n");
          fprintf(fplog, "         H->    \n");
          fprintf(fplog, "         6 10    \n");
          fprintf(fplog, "         H<-    \n");
          fprintf(fplog, "                 \n");
          fprintf(fplog, "         O->    \n");
          fprintf(fplog, "         8 12 \n");
          fprintf(fplog, "         O<-    \n");
          fprintf(fplog, "                 \n"); 
          plumed_error("PluMeD dead with errors: check log file");
  }

  colvar.r_0[count]	 = (real) r_0;
  colvar.type_s[count]   = 7;

  fprintf(fplog, "\n%1i-HBONDS: NN %i MM %i R_0 %lf ", count+1, colvar.nn[count], colvar.mm[count], colvar.r_0[count]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON");
  else                       fprintf(fplog, " PBC OFF");
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  // ##### Type 5 is documented
  fprintf(fplog, "|--HBONDS: WILL COUNT INTRA PROTEIN HBONDS, TYPE = %i (0 = ALL, 1 = BETA ODD, 2 = ALPHA, 3 = BETA EVEN, 4 = NATIVE, 5 = ALL NOT WITHIN 5 ATOMS, 6 = BETA ALL)\n", colvar.type[count]);

  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");  

  snew(colvar.myder[count], colvar.natoms[count]);

  return colvar.natoms[count];
}
