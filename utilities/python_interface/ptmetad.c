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
#include <assert.h>

// This routine has to be executed on the slave nodes to help ptmetad_vbias.
// The call to non-local routines (gmx_sum, hills_engine) HAS TO BE IN THE SAME ORDER as in ptmetad_vbias.
void ptmetad_helper(){
  real s[nconst_max],s2[nconst_max];
  int icv,ncv;
#if defined(PLUMED_GROMACS4) || defined(PLUMED_GROMACS45)
// only slaves should enter this routine
  assert(!mtd_data.ionode);

  ncv=colvar.nconst;

  if(PAR(mtd_data.mcr) && logical.parallel_hills) {
    for(icv=0;icv<ncv;icv++) s[icv]=0.0;
    for(icv=0;icv<ncv;icv++) s2[icv]=0.0;
    plumed_sum(&mtd_data,ncv,s);
    plumed_sum(&mtd_data,ncv,s2);
// returned values are discarded, since they are relevant only for the master node (which calculates the acceptance)
    hills_engine(s,NULL);
    hills_engine(s2,NULL);
  }
#endif
};

void ptmetad_vbias(int newrepl,real*bias,real*biasx){
#if defined(PLUMED_GROMACS) || defined(OPEP) 
// only master process enters here; slaves go to ptmetad_helper
// ATTENTION: non-local calls (hills_engine and gmx_sum) has to be in the same order
  real local_bias; // bias with my hamiltonian in my replica;
  real local_biasx; // bias with my hamiltonian in replica newrepl;
  int icv,irep;
  int ncv,nrep;
  int myrep;
  real **all_ss0;

  ncv=colvar.nconst;
  nrep=mtd_data.nrepl;
  myrep=mtd_data.repl;

// if no exchange has to be done (newrepl<0), we set it to myrep.
// in this way, the bias is computed even when 
  if(newrepl<0) newrepl=myrep;

// array with positions of all replicas
  all_ss0=float_2d_array_alloc(nrep,ncv);

// a few assertions for debugging:
  assert(newrepl<nrep); // inrange
  assert(myrep>=0); assert(myrep<nrep); // in range
  assert(colvar.ss0); // allocated
  assert(all_ss0); // allocated

// sharing current CVs value
  for(icv=0;icv<ncv;icv++) for(irep=0;irep<nrep;irep++) all_ss0[irep][icv] = 0.;
  for(icv=0;icv<ncv;icv++)                              all_ss0[myrep][icv] = colvar.ss0[icv];
  plumed_intersum(&mtd_data,nrep*ncv,& all_ss0[0][0]);

#if defined(PLUMED_GROMACS4) || defined(PLUMED_GROMACS45)
// sharing the cv values with the slaves in ptmetad_helper()
   if(logical.parallel_hills) {
     plumed_sum(&mtd_data,ncv,colvar.ss0);
     plumed_sum(&mtd_data,ncv,all_ss0[newrepl]);
   }
#endif

// calculating the bias, walls and external potential
  assert(all_ss0[newrepl]);
  local_bias=hills_engine(colvar.ss0,NULL)+soft_walls_engine(colvar.ss0,NULL)+ext_forces_engine(colvar.ss0,&extpot,NULL)
             +steer_engine(colvar.ss0,NULL);
  local_biasx=hills_engine(all_ss0[newrepl],NULL)+soft_walls_engine(all_ss0[newrepl],NULL)+ext_forces_engine(all_ss0[newrepl],&extpot,NULL)
             +steer_engine(all_ss0[newrepl],NULL);

// we don't need it anymore
  free_2dr_array_alloc(all_ss0,nrep);
  
// sharing the bias
  for(irep=0;irep<nrep;irep++) bias[irep]=0.;
  for(irep=0;irep<nrep;irep++) biasx[irep]=0.;
  bias[myrep]=local_bias;
  biasx[myrep]=local_biasx;
  plumed_intersum(&mtd_data,nrep, & bias[0]);
  plumed_intersum(&mtd_data,nrep, & biasx[0]);
#endif
};

// In this implementation, all the replicas have to call this routine;
void ptmetad_exchfluct(int newrepl)
{
#if defined(PLUMED_GROMACS) || defined(OPEP) 
  int i,j;
  real **all_Mss0;
  real **all_M2ss0;
  all_Mss0=float_2d_array_alloc(mtd_data.nrepl,colvar.nconst);
  all_M2ss0=float_2d_array_alloc(mtd_data.nrepl,colvar.nconst);
  for(i=0;i<colvar.nconst;i++) for(j=0;j<mtd_data.nrepl;j++) all_Mss0[j][i] = 0.;
  for(i=0;i<colvar.nconst;i++) for(j=0;j<mtd_data.nrepl;j++) all_M2ss0[j][i] = 0.;
  for(i=0;i<colvar.nconst;i++) all_Mss0[i][mtd_data.repl] = colvar.Mss0[i];
  for(i=0;i<colvar.nconst;i++) all_M2ss0[i][mtd_data.repl] = colvar.M2ss0[i];
  plumed_intersum(&mtd_data,mtd_data.nrepl*colvar.nconst, & all_M2ss0[0][0]);
  if(newrepl>=0){
    for(i=0;i<colvar.nconst;i++) colvar.Mss0[i]=all_Mss0[i][newrepl];
    for(i=0;i<colvar.nconst;i++) colvar.M2ss0[i]=all_M2ss0[i][newrepl];
  };
  free_2dr_array_alloc(all_Mss0,mtd_data.nrepl);
  free_2dr_array_alloc(all_M2ss0,mtd_data.nrepl);
#endif
}

