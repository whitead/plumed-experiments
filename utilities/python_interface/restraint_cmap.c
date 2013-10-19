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

void  PREFIX cmap_restraint(int i_c, struct mtd_data_s *mtd_data) {

  int    iat, i, tot_con;
  struct sz_data *pmy_sz;
  struct cmap_inpack inpack;
  real   cmap;
  int    jj,k,j,ii,jjj,iii;
  int    tot;
  real   tmp4_r0_0,tmp4_r0_1,tmp4_r0_2;
  real   tmp2_r0,tmp3_r0, dist_r0;
  real   tmp1,tmp2,tmp3;
  real   pow_P, pow_Q, R01;
  int    P1,Q1;
  rvec   rij;


  cmap = 0.;
  pmy_sz=&my_sz_list[ic_to_sz[i_c]];
                                                         
  for (i=0;i<colvar.natoms[i_c];i++){
               iat = colvar.cvatoms[i_c][i];
               inpack.r0[i][0] = mtd_data->pos[iat][0];
               inpack.r0[i][1] = mtd_data->pos[iat][1];
               inpack.r0[i][2] = mtd_data->pos[iat][2];
  }

  cmap_running(i_c, &inpack,&pmy_sz->my_cmap_pack);

  tot_con=pmy_sz->my_cmap_pack.number+pmy_sz->my_cmap_pack.gnumber;        
  for(i=0;i<tot_con;i++) cmap += inpack.cmap[i];

  colvar.ss0[i_c] = cmap; 

// derivatives calculation
  for(i=0;i<colvar.natoms[i_c];i++) {
     colvar.myder[i_c][i][0] = 0.;
     colvar.myder[i_c][i][1] = 0.; 
     colvar.myder[i_c][i][2] = 0.;
  }

  for(j=0;j<pmy_sz->my_cmap_pack.number;j++){
   if(pmy_sz->my_cmap_pack.weight[j]!=0){
     ii=pmy_sz->my_cmap_pack.index_from1[j];
     jj=pmy_sz->my_cmap_pack.index_from2[j];

     if(colvar.cell_pbc[i_c]){
       minimal_image(inpack.r0[ii], inpack.r0[jj], &dist_r0, rij);
     } else {
       rij[0] = inpack.r0[ii][0]-inpack.r0[jj][0];
       rij[1] = inpack.r0[ii][1]-inpack.r0[jj][1];
       rij[2] = inpack.r0[ii][2]-inpack.r0[jj][2];
       dist_r0=sqrt(pow2(inpack.r0[ii][0]-inpack.r0[jj][0])+
                 pow2(inpack.r0[ii][1]-inpack.r0[jj][1])+
                 pow2(inpack.r0[ii][2]-inpack.r0[jj][2]));
     };

     R01=pmy_sz->my_cmap_pack.r0[j];
     P1=pmy_sz->my_cmap_pack.nn[j];
     Q1=pmy_sz->my_cmap_pack.nd[j];

     if(fabs(dist_r0/R01-1.0)<0.00001){
           colvar.myder[i_c][ii][0]+=rij[0]*P1*(P1-Q1)/Q1;
           colvar.myder[i_c][ii][1]+=rij[1]*P1*(P1-Q1)/Q1;
           colvar.myder[i_c][ii][2]+=rij[2]*P1*(P1-Q1)/Q1;
           colvar.myder[i_c][jj][0]-=rij[0]*P1*(P1-Q1)/Q1;
           colvar.myder[i_c][jj][1]-=rij[1]*P1*(P1-Q1)/Q1;
           colvar.myder[i_c][jj][2]-=rij[2]*P1*(P1-Q1)/Q1;
     } else {
           power(dist_r0/R01,P1,Q1,&pow_P,&pow_Q);

           tmp2_r0=(Q1*pow_Q*(1.-pow_P)-P1*pow_P*(1.-pow_Q))*R01/dist_r0*pmy_sz->my_cmap_pack.weight[j];
           tmp3_r0=R01*(1.-pow_Q)*(1.-pow_Q);

           tmp4_r0_0=rij[0]/dist_r0;
           tmp4_r0_1=rij[1]/dist_r0;
           tmp4_r0_2=rij[2]/dist_r0;

           tmp1=tmp2_r0*tmp4_r0_0/tmp3_r0;
           tmp2=tmp2_r0*tmp4_r0_1/tmp3_r0;
           tmp3=tmp2_r0*tmp4_r0_2/tmp3_r0;
           colvar.myder[i_c][ii][0]+=tmp1;
           colvar.myder[i_c][ii][1]+=tmp2;
           colvar.myder[i_c][ii][2]+=tmp3;
           colvar.myder[i_c][jj][0]-=tmp1;
           colvar.myder[i_c][jj][1]-=tmp2;
           colvar.myder[i_c][jj][2]-=tmp3;
     }
   }
  }

// group contacts

  if(pmy_sz->my_cmap_pack.logical_group){
   for(j=0;j<pmy_sz->my_cmap_pack.gnumber;j++){

    i=j+pmy_sz->my_cmap_pack.number;
    ii=pmy_sz->my_cmap_pack.index1[i];
    jj=pmy_sz->my_cmap_pack.index2[i];

    if(colvar.cell_pbc[i_c]){
      minimal_image(pmy_sz->my_cmap_pack.group.rcm[ii], pmy_sz->my_cmap_pack.group.rcm[jj], &dist_r0, rij);
    } else {
      rij[0] = pmy_sz->my_cmap_pack.group.rcm[ii][0]-pmy_sz->my_cmap_pack.group.rcm[jj][0];
      rij[1] = pmy_sz->my_cmap_pack.group.rcm[ii][1]-pmy_sz->my_cmap_pack.group.rcm[jj][1];
      rij[2] = pmy_sz->my_cmap_pack.group.rcm[ii][2]-pmy_sz->my_cmap_pack.group.rcm[jj][2];
      dist_r0=sqrt(pow2(pmy_sz->my_cmap_pack.group.rcm[ii][0]-pmy_sz->my_cmap_pack.group.rcm[jj][0])+
                   pow2(pmy_sz->my_cmap_pack.group.rcm[ii][1]-pmy_sz->my_cmap_pack.group.rcm[jj][1])+
                   pow2(pmy_sz->my_cmap_pack.group.rcm[ii][2]-pmy_sz->my_cmap_pack.group.rcm[jj][2]));
    };


   R01=pmy_sz->my_cmap_pack.r0[i];
   P1=pmy_sz->my_cmap_pack.nn[i];
   Q1=pmy_sz->my_cmap_pack.nd[i];

   if(fabs(dist_r0/R01-1.0)<0.00001){
    for(jjj=0;jjj<pmy_sz->my_cmap_pack.group.numatom[ii];jjj++){
     iii=pmy_sz->my_cmap_pack.group.index_to_list[ii][jjj];
     colvar.myder[i_c][iii][0]+=rij[0]*P1*(P1-Q1)/Q1/((real) pmy_sz->my_cmap_pack.group.numatom[ii]);
     colvar.myder[i_c][iii][1]+=rij[1]*P1*(P1-Q1)/Q1/((real) pmy_sz->my_cmap_pack.group.numatom[ii]);
     colvar.myder[i_c][iii][2]+=rij[2]*P1*(P1-Q1)/Q1/((real) pmy_sz->my_cmap_pack.group.numatom[ii]);
    }
    for(jjj=0;jjj<pmy_sz->my_cmap_pack.group.numatom[jj];jjj++){
     iii=pmy_sz->my_cmap_pack.group.index_to_list[jj][jjj];
     colvar.myder[i_c][iii][0]-=rij[0]*P1*(P1-Q1)/Q1/((real) pmy_sz->my_cmap_pack.group.numatom[jj]);
     colvar.myder[i_c][iii][1]-=rij[1]*P1*(P1-Q1)/Q1/((real) pmy_sz->my_cmap_pack.group.numatom[jj]);
     colvar.myder[i_c][iii][2]-=rij[2]*P1*(P1-Q1)/Q1/((real) pmy_sz->my_cmap_pack.group.numatom[jj]);
    }

   }else{

    power(dist_r0/R01,P1,Q1,&pow_P,&pow_Q);

    tmp2_r0=(Q1*pow_Q*(1.-pow_P)-P1*pow_P*(1.-pow_Q))*R01/dist_r0*pmy_sz->my_cmap_pack.weight[i];
    tmp3_r0=R01*(1.-pow_Q)*(1.-pow_Q);

    tmp4_r0_0=rij[0]/dist_r0;
    tmp4_r0_1=rij[1]/dist_r0;
    tmp4_r0_2=rij[2]/dist_r0;

    tmp1=tmp2_r0*tmp4_r0_0/tmp3_r0;
    tmp2=tmp2_r0*tmp4_r0_1/tmp3_r0;
    tmp3=tmp2_r0*tmp4_r0_2/tmp3_r0;

    for(jjj=0;jjj<pmy_sz->my_cmap_pack.group.numatom[ii];jjj++){
     iii=pmy_sz->my_cmap_pack.group.index_to_list[ii][jjj];
     colvar.myder[i_c][iii][0]+=tmp1/((real) pmy_sz->my_cmap_pack.group.numatom[ii]);
     colvar.myder[i_c][iii][1]+=tmp2/((real) pmy_sz->my_cmap_pack.group.numatom[ii]);
     colvar.myder[i_c][iii][2]+=tmp3/((real) pmy_sz->my_cmap_pack.group.numatom[ii]);
    }
    for(jjj=0;jjj<pmy_sz->my_cmap_pack.group.numatom[jj];jjj++){
     iii=pmy_sz->my_cmap_pack.group.index_to_list[jj][jjj];
     colvar.myder[i_c][iii][0]-=tmp1/((real) pmy_sz->my_cmap_pack.group.numatom[jj]);
     colvar.myder[i_c][iii][1]-=tmp2/((real) pmy_sz->my_cmap_pack.group.numatom[jj]);
     colvar.myder[i_c][iii][2]-=tmp3/((real) pmy_sz->my_cmap_pack.group.numatom[jj]);
    }
   }
  }
 }


  return;
}

// ------------------------------------------------------------------------------------------------

int PREFIX read_cmap(char **word, int count,t_plumed_input *input,FILE *fplog)
{

  int i,iw;
  double sigma = 0.0;
  char file_maps[129];
  char file_group[129];
  struct sz_data *my_sz;
  int help;

  help=0; 
  colvar.cell_pbc[count]=1; // default is PBC

  my_sz=&(my_sz_list[nsz]);
  ic_to_sz[count]=nsz;

  my_sz->my_cmap_pack.logical_group = 0; // no GROUP

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &sigma);
             colvar.delta_r[count]  = (real) sigma; }
  iw = seek_word(word,"INDEX");
  if(iw>=0) { sscanf(word[iw+1],"%s", file_maps);
  } else {
    fprintf(fplog,"|- NEEDED INDEX KEYWORD FOR CMAP\n");
    help=1;
  }
  iw=seek_word(word,"GROUP");
  if(iw>=0) {
   sscanf(word[iw+1],"%s",file_group);
   my_sz->my_cmap_pack.logical_group = 1;
  }
  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}

  if(help){
         fprintf(fplog,"|- CMAP SYNTAX:\n");
         fprintf(fplog,"|- INDEX              : contact definition file \n");
         fprintf(fplog,"|- GROUP              : group definition file \n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"|- e.g.\n");
         fprintf(fplog,"|- \n");
         fprintf(fplog,"|- CMAP INDEX cmap.index { GROUP cmap.group } { NOPBC } SIGMA 1.0 \n");
         fprintf(fplog,"|- \n");
         plumed_error("PluMeD dead with errors: check log file");
  } 

  
  fprintf(fplog, "\n%1i-CMAP: INDEX file %s ", count+1, file_maps);
  if(my_sz->my_cmap_pack.logical_group) fprintf(fplog, " GROUP file %s ", file_group); 
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON ");
  else                       fprintf(fplog, " PBC OFF ");
  if (logical.do_hills){
        if (colvar.delta_r[count]>0){
                 fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
        }
  }
  else fprintf(fplog,"\n");


  read_sz_map(my_sz,file_maps,file_maps,file_group,0,fplog);

  colvar.natoms[count]   = my_sz->my_cmap_pack.atoms; 

  snew(colvar.myder[count], colvar.natoms[count]);
  snew(colvar.cvatoms[count], colvar.natoms[count]);

  for(i=0;i<colvar.natoms[count];i++)
    colvar.cvatoms[count][i] = my_sz->my_cmap_pack.list[i];
  
  nsz++;

  if(nsz==NMAX_PATH) plumed_error("TOO MANY PATH CVS. Increase NMAX_PATH in metadyn.h and recompile");

  fprintf(fplog,"\n");

  return colvar.natoms[count];
}
