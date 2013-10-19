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

void PREFIX alpharmsd_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int firstAtom, secondAtom, i, j, ix;
  int nn = colvar.nn[i_c], mm = colvar.mm[i_c];
  rvec rij;
  real num, iden, mod_rij, rdist, func, dfunc, rNdist, rMdist;
  real ncoord, r_0 = colvar.r_0[i_c], d_0 = colvar.d_0[i_c];
  real threshold;
  threshold=pow(0.00001,1./(nn-mm));

  int ires, nres, iat, add_grad;
  int atom_list[30];
  real dist, g[30][3], rmsd;

  nres=colvar.natoms[i_c]/5;
  ncoord=0.;
  for (i=0;i<colvar.natoms[i_c];i++) { for (ix=0;ix<3;ix++) { colvar.myder[i_c][i][ix]=0.; } }

  // ires:   0 1 ... = (N CA C O CB) (N CA C O CB) ...
  // iat:    0 1 2 3 4 5 6 ... = N CA C O CB N CA ...
  // loop over all segments of 6 residues in the protein chain
  for (ires=0;ires<nres-5;ires++) { 
    // copying atom indexes from 6 residues (6*5 atoms) in atom_list
    i=0;
    for (iat=ires*5;iat<(ires+6)*5;iat++) { 
      atom_list[i]=iat;
//      if(colvar.it==0){printf("ires1 = %d atom_list %d = %d\n",ires,i,colvar.cvatoms[i_c][atom_list[i]]);} // DEBUG
      i++;
    }
    // reset
    rmsd=0.;
    for (i=0;i<30;i++) { for (ix=0;ix<3;ix++) { g[i][ix]=0.; } }
    // rmsd between distance matrices
    for (i=0;i<29;i++) {
      firstAtom=colvar.cvatoms[i_c][atom_list[i]];
      for (j=i+1;j<30;j++) {
        secondAtom=colvar.cvatoms[i_c][atom_list[j]];
        // skip covalent bonds
        if (ref_dist_mat.alpha[i][j]<=0.) continue;
        if(colvar.cell_pbc[i_c]){
          minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
        } else {
          rij[0] = mtd_data->pos[firstAtom][0]-mtd_data->pos[secondAtom][0];
          rij[1] = mtd_data->pos[firstAtom][1]-mtd_data->pos[secondAtom][1];
          rij[2] = mtd_data->pos[firstAtom][2]-mtd_data->pos[secondAtom][2];
          mod_rij = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
        }
        dist = mod_rij-ref_dist_mat.alpha[i][j];
//        if(colvar.it==0&&dist>0.2){printf("i j ai aj rij rij0 = %d %d %d %d %8.3f %8.3f\n",i,j,colvar.cvatoms[i_c][atom_list[i]],colvar.cvatoms[i_c][atom_list[j]],mod_rij,ref_dist_mat.alpha[i][j]);} // DEBUG
        rmsd += dist*dist;
        // store the atom-dependent part of the gradients
        for (ix=0;ix<3;ix++) {
          g[i][ix]+=dist*rij[ix]/mod_rij;
          g[j][ix]-=dist*rij[ix]/mod_rij;
        }
      }
    }
    // normalize by the number of off-diagonal elements in the 30*30 matrix 
    rmsd=sqrt(rmsd/ref_dist_mat.alpha_pairs); 
//    printf("step = %d  ires = %d  RMSD = %8.3f\n",colvar.it,ires,rmsd); // DEBUG
    // switching function
    rdist = (rmsd-d_0)/r_0;
    add_grad=0;
    /* analitic limit of the switching function */
    if(rdist<=0.){
      ncoord+=1.;
      dfunc=0.;
    }else if(rdist>0.999999 && rdist<1.000001){
      ncoord+=nn/mm;
      dfunc=0.5*nn*(nn-mm)/mm;
      add_grad=1;
    }else if(rdist>threshold){
      dfunc=0.;
    }else{
      rNdist = pow(rdist, nn-1);
      rMdist = pow(rdist, mm-1);
      num = 1.-rNdist*rdist;
      iden = 1./(1.-rMdist*rdist);
      func = num*iden;
      ncoord += func;
      dfunc = (-nn*rNdist*iden)+(func*(iden*mm)*rMdist);
      add_grad=1;
    }
    if (add_grad==1) {
      dfunc=dfunc/(r_0*rmsd*ref_dist_mat.alpha_pairs);
      // total gradient
      for (i=0;i<30;i++) {
        for(ix=0;ix<3;ix++) {
          colvar.myder[i_c][atom_list[i]][ix] += dfunc*g[i][ix];
        }
      }
    }
  } // loop on ires
  colvar.ss0[i_c] = ncoord;

}

//--------------------------------------------------------------------------------------------------------------

int PREFIX read_alpharmsd(char **word, int count, t_plumed_input *input, FILE *fplog)
{
// Coordinates N CA C O CB (Angstrom) of alpha 6-res block from the representative pdbs of
// each of the 5 architecture entries in class "mainly alpha" of CATH database
// (1oai 1mz9 1jdh 1ppr 1h12)
  static rvec ref_alpha[30]={   
    { 0.733,  0.519,  5.298}, // N    i
    { 1.763,  0.810,  4.301}, // CA
    { 1.527, -0.045,  3.053}, // C
    { 1.646,  0.436,  1.928}, // O
    { 3.166,  0.543,  4.881}, // CB
    { 1.180, -1.312,  3.254}, // N    i+1
    { 0.924, -2.203,  2.126}, // CA
    {-0.239, -1.711,  1.261}, // C
    {-0.190, -1.815,  0.032}, // O
    { 0.650, -3.626,  2.626}, // CB
    {-1.280, -1.172,  1.891}, // N    i+2
    {-2.416, -0.661,  1.127}, // CA
    {-1.964,  0.529,  0.276}, // C
    {-2.364,  0.659, -0.880}, // O
    {-3.548, -0.217,  2.056}, // CB
    {-1.130,  1.391,  0.856}, // N    i+3
    {-0.620,  2.565,  0.148}, // CA
    { 0.231,  2.129, -1.032}, // C
    { 0.179,  2.733, -2.099}, // O
    { 0.228,  3.439,  1.077}, // CB
    { 1.028,  1.084, -0.833}, // N    i+4
    { 1.872,  0.593, -1.919}, // CA
    { 1.020,  0.020, -3.049}, // C
    { 1.317,  0.227, -4.224}, // O
    { 2.850, -0.462, -1.397}, // CB
    {-0.051, -0.684, -2.696}, // N    i+5
    {-0.927, -1.261, -3.713}, // CA
    {-1.663, -0.171, -4.475}, // C
    {-1.916, -0.296, -5.673}, // O
    {-1.933, -2.219, -3.074}  // CB
  };

  int i, iw, iat, j, help, ix;
  double r_0, d_0, d;
  double delta = 0.0;
  double angstrom_scale;
  char string[400];
  real threshold, value;

  help=0;
  d_0=0.;

  colvar.cell_pbc[count]=1; // default is PBC

  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR ALPHARMSD\n"); help=1;}

  if(colvar.natoms[count]%5!=0){ fprintf(fplog,"|- ERROR: TOTAL NUMBER OF ATOMS MUST BE MULTIPLE OF 5\n"); help=1; }

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }
  iw=seek_word(word,"NN");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.nn[count]); } else { fprintf(fplog,"|- NEEDED NN KEYWORD FOR ALPHARMSD\n"); help=1;}
  iw=seek_word(word,"MM");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.mm[count]);} else { fprintf(fplog,"|- NEEDED MM KEYWORD FOR ALPHARMSD\n"); help=1;}
  iw=seek_word(word,"R_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &r_0); } else { fprintf(fplog,"|- NEEDED R_0 KEYWORD FOR ALPHARMSD\n"); help=1;}
  iw=seek_word(word,"D_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &d_0); }
  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}
  iw=seek_word(word,"ANGSTROM_SCALE");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &angstrom_scale); } else { fprintf(fplog,"|- NEEDED ANGSTROM_SCALE KEYWORD FOR ALPHARMSD\n"); help=1;}

  if(help){
          fprintf(fplog, "\n-ALPHARMSD CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "ALPHARMSD LIST <g1> NN 8 MM 12 R_0 0.08 SIGMA 1.0 ANGSTROM_SCALE 0.1 NOPBC\n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "            7     9    26    27    11 \n");
          fprintf(fplog, "           28    30    45    46    32 \n");
          fprintf(fplog, "           47    49    56    57    51 \n");
          fprintf(fplog, "           58    60    71    72    62 \n");
          fprintf(fplog, "           73    75    88    89    77 \n");
          fprintf(fplog, "           90    92   100   101    94 \n");
          fprintf(fplog, "          102   104   120   121   106 \n");
          fprintf(fplog, "          122   124   136   137   126 \n");
          fprintf(fplog, "          138   140   147   148   142 \n");
          fprintf(fplog, "          149   151   163   164   153 \n");
          fprintf(fplog, "          165   167   183   184   169 \n");
          fprintf(fplog, "          185   187   190   191   189 \n");
          fprintf(fplog, "          192   194   209   210   196    \n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "       \n");
          fprintf(fplog, " (where each row contains N CA C O CB of consecutive residues)\n");
          fprintf(fplog, "       \n");
          plumed_error("PluMed dead with errors: check log file");
  }

  colvar.r_0[count]      = (real) r_0;
  colvar.d_0[count]      = (real) d_0;
  colvar.type_s[count]   = 37;

  fprintf(fplog, "%1i-ALPHARMSD; ATOMS INVOLVED: %i; ", count+1, colvar.natoms[count]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON");
  else                       fprintf(fplog, " PBC OFF");
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  fprintf(fplog, "|--FUNCTIONAL FORM: (1-((dist_mat_rmsd-d_0)/r_0)^n) / (1-((dist_mat_rmsd-d_0)/r_0)^m) \n");
  fprintf(fplog, "|--PARAMETERS: n= %i m= %i r_0= %f d_0= %f\n", colvar.nn[count], colvar.mm[count], colvar.r_0[count], colvar.d_0[count]);
  threshold=pow(0.00001,1./(colvar.nn[count]-colvar.mm[count]));
  value=(1.-pow(threshold,colvar.nn[count]))/(1.-pow(threshold,colvar.mm[count]));
  fprintf(fplog, "|--CUTOFF VALUE: %f\n",value);
  fprintf(fplog, "|--CUTOFF DISTANCE: %f\n",threshold*r_0+d_0);

  iat=0;
  fprintf(fplog,"|- SET MEMBERS: ");
  for(i=0;i<colvar.natoms[count];i++){fprintf(fplog," %d ",colvar.cvatoms[count][i]+1);if((i+1)%20==0)fprintf(fplog,"\n               ");}fprintf(fplog,"\n\n");

  fprintf(fplog,"|- IMPORTANT NOTE: THIS CV REQUIRES ALIGN_ATOMS (FOR ALL BACKBONE) ON SOME CODES LIKE GROMACS4 !!!\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  // storing reference distance matrix (converted from angstrom to angstrom*angstrom_scale)
  ref_dist_mat.alpha_pairs=0.; // pairs without covalent bond
  for (i=0;i<29;i++) {
    for (j=i+1;j<30;j++) {
      d=0.;
      for (ix=0;ix<3;ix++) {
        d += pow((ref_alpha[i][ix]-ref_alpha[j][ix]),2);
      }
      d=sqrt(d);
      // here d is in angstrom: covalent bonds among N CA C O CB are maximum 1.54 A
      if(d>1.7) {
        // convert d in length units of the MD program
        d=d*angstrom_scale;
        ref_dist_mat.alpha_pairs+=1.;
      } else {
        // exclude covalent bonds
        d=-999.;
      }
      ref_dist_mat.alpha[i][j]=d;
      ref_dist_mat.alpha[j][i]=d;
    }
  }
 
  return colvar.natoms[count]; 
}
