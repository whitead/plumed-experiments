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

void PREFIX parabetarmsd_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int firstAtom, secondAtom, i, j, ix;
  int nn = colvar.nn[i_c], mm = colvar.mm[i_c];
  rvec rij;
  real num, iden, mod_rij, rdist, func, dfunc, rNdist, rMdist;
  real ncoord, r_0 = colvar.r_0[i_c], d_0 = colvar.d_0[i_c];
  real threshold;
  threshold=pow(0.00001,1./(nn-mm));

  int ires, jres, nres, iat, add_grad, iblock;
  int atom_list[30];
  real dist_mat[30][30], dist_mat_rij[30][30][3];
  real dist, g[30][3], rmsd, cutoff;
  cutoff=colvar.realpar[i_c][0][0];

  nres=colvar.natoms[i_c]/5;
  ncoord=0.;
  for (i=0;i<colvar.natoms[i_c];i++) { for (ix=0;ix<3;ix++) { colvar.myder[i_c][i][ix]=0.; } }

  // ires:   0 1 ... = (N CA C O CB) (N CA C O CB) ...
  // iat:    0 1 2 3 4 5 6 ... = N CA C O CB N CA ...
  // loop over residue triplets {i,j}{i+1,j+1}{i+2,j+2} 
  //  with (j-i)>5 (this leaves min 3 residues in the connecting chain)
  for (ires=0;ires<nres-8;ires++) { 
  for (jres=ires+6;jres<nres-2;jres++) { 
    // copying atom indexes from 3 pairs of residues (6*5 atoms) in atom_list
    // order of H-bonded residues in ref_parabeta[]: 1-4 2-5 3-6
    for(i=0;i<15;i++) {
      atom_list[i]=ires*5+i;
      atom_list[i+15]=jres*5+i;
    }
//    for(i=0;i<30;i++){printf("ires1,jres1 = %d %d atom_list %d = %d\n",ires,jres,i,colvar.cvatoms[i_c][atom_list[i]]);} // DEBUG
    // check if it's worth to compute the rmsd (if segments are not too far)
    firstAtom=colvar.cvatoms[i_c][atom_list[6]];   // CA at the center of first segment
    secondAtom=colvar.cvatoms[i_c][atom_list[21]]; // CA at the center of second segment
    if(colvar.cell_pbc[i_c]){
      minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
    } else {
      rij[0] = mtd_data->pos[firstAtom][0]-mtd_data->pos[secondAtom][0];
      rij[1] = mtd_data->pos[firstAtom][1]-mtd_data->pos[secondAtom][1];
      rij[2] = mtd_data->pos[firstAtom][2]-mtd_data->pos[secondAtom][2];
      mod_rij = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
    }
    if (mod_rij>cutoff) continue; // this can save a lot of computer time!
    // rmsd between distance matrices
    for (i=0;i<29;i++) {
      firstAtom=colvar.cvatoms[i_c][atom_list[i]];
      for (j=i+1;j<30;j++) {
        secondAtom=colvar.cvatoms[i_c][atom_list[j]];
        // skip covalent bonds
        if (ref_dist_mat.parabeta[0][i][j]<=0.||ref_dist_mat.parabeta[1][i][j]<=0.) continue;
        if(colvar.cell_pbc[i_c]){
          minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
        } else {
          rij[0] = mtd_data->pos[firstAtom][0]-mtd_data->pos[secondAtom][0];
          rij[1] = mtd_data->pos[firstAtom][1]-mtd_data->pos[secondAtom][1];
          rij[2] = mtd_data->pos[firstAtom][2]-mtd_data->pos[secondAtom][2];
          mod_rij = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
        }
        dist_mat[i][j]=mod_rij;
        for (ix=0;ix<3;ix++) { dist_mat_rij[i][j][ix]=rij[ix]; }
      }
    }
    // loop over the two types of parabeta blocks
    for (iblock=0;iblock<2;iblock++) {
      // reset
      rmsd=0.;
      for (i=0;i<30;i++) { for (ix=0;ix<3;ix++) { g[i][ix]=0.; } }
      for (i=0;i<29;i++) {
        for (j=i+1;j<30;j++) {
          if (ref_dist_mat.parabeta[0][i][j]<=0.||ref_dist_mat.parabeta[1][i][j]<=0.) continue;
          dist = dist_mat[i][j]-ref_dist_mat.parabeta[iblock][i][j];
//          printf("i j ai aj rij rij0 = %d %d %d %d %8.3f %8.3f\n",i,j,colvar.cvatoms[i_c][atom_list[i]],colvar.cvatoms[i_c][atom_list[j]],dist_mat[i][j],ref_dist_mat.parabeta[iblock][i][j]); // DEBUG
          rmsd += dist*dist;
          // store the atom-dependent part of the gradients
          for (ix=0;ix<3;ix++) {
            g[i][ix]+=dist*dist_mat_rij[i][j][ix]/dist_mat[i][j];
            g[j][ix]-=dist*dist_mat_rij[i][j][ix]/dist_mat[i][j];
          }
        }
      }
      // normalize by the number of off-diagonal elements in the 30*30 matrix 
      rmsd=sqrt(rmsd/ref_dist_mat.parabeta_pairs[iblock]); 
//      printf("step = %d  ires,jres = %d %d  RMSD = %8.3f\n",colvar.it,ires,jres,rmsd); // DEBUG
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
        dfunc=dfunc/(r_0*rmsd*ref_dist_mat.parabeta_pairs[iblock]);
        // total gradient
        for (i=0;i<30;i++) {
          for(ix=0;ix<3;ix++) {
            // the factor 0.5 cures the double-counting in not-too-short sheets
            colvar.myder[i_c][atom_list[i]][ix] += 0.5*dfunc*g[i][ix];
          }
        }
      }
    } // loop on iblock
  } // loop on jres
  } // loop on ires
  // the factor 0.5 cures the double-counting in not-too-short sheets
  colvar.ss0[i_c] = 0.5*ncoord;

}

//--------------------------------------------------------------------------------------------------------------

int PREFIX read_parabetarmsd(char **word, int count, t_plumed_input *input, FILE *fplog)
{
// Coordinates N CA C O CB (Angstrom) of parabeta 3+3 blocks from the representative pdbs of
// each of the 20 architecture entries in class "mainly beta" of CATH database
// (1bds 1gvk 1h8p 1i5p 1itv 1k7i 1m3y 1n7v 1nh2 1qre
// 1rg8 1tl2 1w6s 1ylh 2bbk 2dpf 2hnu 2nwf 3sil 4bcl)
// There are 2 different ideal parabeta blocks!
  static rvec ref_parabeta[60]={
    { 1.244, -4.620, -2.127}, // N    i
    {-0.016, -4.500, -1.395}, // CA
    {-0.287, -3.000, -1.301}, // C
    { 0.550, -2.245, -0.822}, // O
    { 0.105, -5.089,  0.024}, // CB
    {-1.445, -2.551, -1.779}, // N    i+1
    {-1.752, -1.130, -1.677}, // CA
    {-2.906, -0.961, -0.689}, // C
    {-3.867, -1.738, -0.695}, // O
    {-2.113, -0.550, -3.059}, // CB
    {-2.774,  0.034,  0.190}, // N    i+2
    {-3.788,  0.331,  1.201}, // CA
    {-4.294,  1.743,  0.937}, // C
    {-3.503,  2.671,  0.821}, // O
    {-3.188,  0.300,  2.624}, // CB
    { 4.746, -2.363,  0.188}, // N    j
    { 3.427, -1.839,  0.545}, // CA
    { 3.346, -0.365,  0.181}, // C
    { 4.237,  0.412,  0.521}, // O
    { 3.135, -1.958,  2.074}, // CB
    { 2.261,  0.013, -0.487}, // N    j+1
    { 2.024,  1.401, -0.875}, // CA
    { 0.914,  1.902,  0.044}, // C
    {-0.173,  1.330,  0.052}, // O
    { 1.489,  1.514, -2.313}, // CB
    { 1.202,  2.940,  0.828}, // N    j+2
    { 0.190,  3.507,  1.718}, // CA
    {-0.229,  4.791,  1.038}, // C
    { 0.523,  5.771,  0.996}, // O
    { 0.772,  3.801,  3.104}, // CB
      {-1.439, -5.122, -1.144}, // N    i
      {-0.816, -3.803, -1.013}, // CA
      {-1.928, -2.770, -0.952}, // C
      {-2.991, -2.970, -1.551}, // O
      { 0.099, -3.509, -2.206}, // CB
      {-1.698, -1.687, -0.215}, // N    i+1
      {-2.681, -0.613, -0.143}, // CA
      {-1.984,  0.681, -0.574}, // C
      {-0.807,  0.921, -0.273}, // O
      {-3.323, -0.477,  1.267}, // CB
      {-2.716,  1.492, -1.329}, // N    i+2
      {-2.196,  2.731, -1.883}, // CA
      {-2.989,  3.949, -1.433}, // C
      {-4.214,  3.989, -1.583}, // O
      {-2.263,  2.692, -3.418}, // CB
      { 2.464, -4.352,  2.149}, // N    j
      { 3.078, -3.170,  1.541}, // CA
      { 2.080, -2.021,  1.639}, // C
      { 0.938, -2.178,  1.225}, // O
      { 3.398, -3.415,  0.060}, // CB
      { 2.525, -0.886,  2.183}, // N    j+1
      { 1.692,  0.303,  2.346}, // CA
      { 2.420,  1.410,  1.608}, // C
      { 3.567,  1.733,  1.937}, // O
      { 1.541,  0.665,  3.842}, // CB
      { 1.758,  1.976,  0.600}, // N    j+2
      { 2.373,  2.987, -0.238}, // CA
      { 1.684,  4.331, -0.148}, // C
      { 0.486,  4.430, -0.415}, // O
      { 2.367,  2.527, -1.720}  // CB
  };

  int i, iw, iat, j, help, ix, iblock;
  double r_0, d_0, d;
  double delta = 0.0;
  double angstrom_scale, strands_cutoff;
  char string[400];
  real threshold, value;

  help=0;
  d_0=0.;

  colvar.cell_pbc[count]=1; // default is PBC
  colvar.realpar[count][0][0]=1000000.; // by default, no cutoff on the distance between beta-strands

  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR PARABETARMSD\n"); help=1;}

  if(colvar.natoms[count]%5!=0){ fprintf(fplog,"|- ERROR: TOTAL NUMBER OF ATOMS MUST BE MULTIPLE OF 5\n"); help=1; }

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }
  iw=seek_word(word,"NN");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.nn[count]); } else { fprintf(fplog,"|- NEEDED NN KEYWORD FOR PARABETARMSD\n"); help=1;}
  iw=seek_word(word,"MM");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.mm[count]);} else { fprintf(fplog,"|- NEEDED MM KEYWORD FOR PARABETARMSD\n"); help=1;}
  iw=seek_word(word,"R_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &r_0); } else { fprintf(fplog,"|- NEEDED R_0 KEYWORD FOR PARABETARMSD\n"); help=1;}
  iw=seek_word(word,"D_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &d_0); }
  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}
  iw=seek_word(word,"ANGSTROM_SCALE");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &angstrom_scale); } else { fprintf(fplog,"|- NEEDED ANGSTROM_SCALE KEYWORD FOR PARABETARMSD\n"); help=1;}
  iw=seek_word(word,"STRANDS_CUTOFF");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &strands_cutoff); }

  if(help){
          fprintf(fplog, "\n-PARABETARMSD CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "PARABETARMSD LIST <g1> NN 8 MM 12 R_0 0.08 SIGMA 1.0 ANGSTROM_SCALE 0.1 STRANDS_CUTOFF 1. NOPBC\n");
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
          plumed_error("PluMeD dead with errors: check log file");
  }

  colvar.r_0[count]      = (real) r_0;
  colvar.d_0[count]      = (real) d_0;
  colvar.realpar[count][0][0] = (real) strands_cutoff;
  colvar.type_s[count]   = 39;

  fprintf(fplog, "%1i-PARABETARMSD; ATOMS INVOLVED: %i; ", count+1, colvar.natoms[count]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON;");
  else                       fprintf(fplog, " PBC OFF;");
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  fprintf(fplog, "|--FUNCTIONAL FORM: (1-((dist_mat_rmsd-d_0)/r_0)^n) / (1-((dist_mat_rmsd-d_0)/r_0)^m) \n");
  fprintf(fplog, "|--PARAMETERS: n= %i m= %i r_0= %f d_0= %f\n", colvar.nn[count], colvar.mm[count], colvar.r_0[count], colvar.d_0[count]);
  threshold=pow(0.00001,1./(colvar.nn[count]-colvar.mm[count]));
  value=(1.-pow(threshold,colvar.nn[count]))/(1.-pow(threshold,colvar.mm[count]));
  fprintf(fplog, "|--CUTOFF VALUE: %f\n",value);
  fprintf(fplog, "|--CUTOFF DISTANCE: %f\n",threshold*r_0+d_0);
  if(colvar.realpar[count][0][0]<1000000.) { fprintf(fplog, "|--STRANDS_CUTOFF DISTANCE: %f\n",colvar.realpar[count][0][0]); }
  fprintf(fplog, "|--ANGSTROM_SCALE: %f\n",angstrom_scale);

  iat=0;
  fprintf(fplog,"|- SET MEMBERS: ");
  for(i=0;i<colvar.natoms[count];i++){fprintf(fplog," %d ",colvar.cvatoms[count][i]+1);if((i+1)%20==0)fprintf(fplog,"\n               ");}fprintf(fplog,"\n\n");

  fprintf(fplog,"|- IMPORTANT NOTE: THIS CV REQUIRES ALIGN_ATOMS (FOR ALL BACKBONE) ON SOME CODES LIKE GROMACS4 !!!\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  // storing reference distance matrix (converted from angstrom to angstrom*angstrom_scale)
  ref_dist_mat.parabeta_pairs[0]=0.; // pairs without covalent bond
  ref_dist_mat.parabeta_pairs[1]=0.; // pairs without covalent bond
  for (iblock=0;iblock<2;iblock++) {
    for (i=0;i<29;i++) {
      for (j=i+1;j<30;j++) {
        d=0.;
        for (ix=0;ix<3;ix++) {
          d += pow((ref_parabeta[i+iblock*30][ix]-ref_parabeta[j+iblock*30][ix]),2);
        }
        d=sqrt(d);
        // here d is in angstrom: covalent bonds among N CA C O CB are maximum 1.54 A
        if(d>1.7) {
          // convert d in length units of the MD program
          d=d*angstrom_scale;
          ref_dist_mat.parabeta_pairs[iblock]+=1.;
        } else {
          // exclude covalent bonds
          d=-999.;
        }
        ref_dist_mat.parabeta[iblock][i][j]=d;
        ref_dist_mat.parabeta[iblock][j][i]=d;
      }
    }
  }
 
  return colvar.natoms[count]; 
}
