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

// This routine calculates, from the variaton of a collective variable, the width
// of the gaussian hills along the CV direction.

inline
int int_floor(real number) {
  return (int) number < 0.0 ? -ceil(fabs(number)) : floor(number);                                                                   
}                                                                                                                                    
                                                                         
void PREFIX hills_adapt()
{
  real fluct, step;
  real Mss02, M2ss0;
  int i_c;

  for(i_c=0;i_c<colvar.nconst;i_c++){
   if(colvar.adapt[i_c].on){ 
// Time for recording fluctuations???
    if((logical.not_same_step) && colvar.it%colvar.adapt[i_c].stride==0 && !firstTime ){
      colvar.Mss0[i_c] += colvar.ss0[i_c];
      colvar.M2ss0[i_c] += colvar.ss0[i_c]*colvar.ss0[i_c];
    }
// Time for evaluating width???    
    if((logical.not_same_step) && colvar.it%colvar.adapt[i_c].block==0 && !firstTime ){
      step = (real) colvar.adapt[i_c].stride / colvar.adapt[i_c].block;      
      M2ss0 = colvar.M2ss0[i_c]*step;
      Mss02 = colvar.Mss0[i_c]*colvar.Mss0[i_c]*step*step;  
      if(M2ss0>Mss02) fluct = sqrt(M2ss0-Mss02);
      else fluct = 0.;
      colvar.delta_r[i_c] = colvar.adapt[i_c].widthmultiplier*fluct;
      if(colvar.delta_r[i_c] > colvar.adapt[i_c].sup) colvar.delta_r[i_c] = colvar.adapt[i_c].sup;
      if(colvar.delta_r[i_c] < colvar.adapt[i_c].inf) colvar.delta_r[i_c] = colvar.adapt[i_c].inf;
      colvar.M2ss0[i_c] = 0.;
      colvar.Mss0[i_c] = 0.;
    }
   }
  }
}

//------------------------------------------------------------------------------------------

void PREFIX hills_push(struct mtd_data_s *mtd_data,real ww, real* ss,real* delta)
{
  int icv,ncv;
  static FILE *file=NULL;
  real inv_ss0[nconst_max];
  int wall;
  static int first=1;
  int nactive;
// INVERSION VARS
  int j, i_c;
  real tmps,tmp,tmpd,tmpF1,tmpF2,tmpF3,smooth_func;  
// END INVERSION VARS

  ncv=colvar.nconst;
  if(hills.n_hills+10>hills.ntothills) hills_reallocate(mtd_data); 

  hills.ww[hills.n_hills] = ww;
  for(icv=0;icv<ncv;icv++) hills.ss0_t[hills.n_hills][icv] = ss[icv];	// new hill center
  for(icv=0;icv<ncv;icv++) colvar.delta_s[hills.n_hills][icv] = delta[icv];       // new hill width




// PRINT
  if(!file) file = fopen((mtd_data->ionode?mtd_data->hilfilen:"/dev/null"), "a");
   
// header: list active CVs (useful e.g. for bias-exchange post-processing)
  if(first) {
    if(strlen(colvar.hills_label)>0){ 
      nactive=0;
      for(i_c=0;i_c<colvar.nconst;i_c++){ if(colvar.on[i_c]) nactive++; }
      fprintf(file, "#! ACTIVE %d",nactive);
      if(nactive>0) {
        for(i_c=0;i_c<colvar.nconst;i_c++){ if(colvar.on[i_c]) fprintf(file, " %d",i_c+1); }
      }
      fprintf(file, " %s",colvar.hills_label); 
      fprintf(file,"\n");
    }
    first=0;
  }

  fprintf(file, "%10.3f   ", mtd_data->time);
  for(icv=0;icv<ncv;icv++) if(colvar.on[icv]) fprintf(file, "%14.9f   ", hills.ss0_t[hills.n_hills][icv]);
  for(icv=0;icv<ncv;icv++) if(colvar.on[icv]) fprintf(file, "%14.9f   ", colvar.delta_s[hills.n_hills][icv]);
  if(logical.welltemp){
   fprintf(file, "%14.9f   %4.3f \n", hills.ww[hills.n_hills]*colvar.wfactor/(colvar.wfactor-1.0)/mtd_data->eunit,colvar.wfactor);
  } else {
   fprintf(file, "%14.9f   %4.3f \n", hills.ww[hills.n_hills]/mtd_data->eunit,0.0);
  } 
  
  hills.n_hills++;                              // hills added
// flush all the time standalone: not a big deal... 
#ifdef STANDALONE 
  fclose(file);
#endif
  if(!logical.do_walkers) hills.read=hills.n_hills;

// WALLS
  wall = 0;

  for(icv=0;icv<ncv;icv++) {
    inv_ss0[icv] = colvar.ss0[icv];
    if(logical.ureflect[icv] && colvar.ss0[icv]>(cvw.upper[icv]-colvar.delta_r[icv]) && 
       colvar.ss0[icv]<cvw.upper[icv] && colvar.on[icv]) { inv_ss0[icv] = 2.*cvw.upper[icv]-colvar.ss0[icv]; wall=1;}
    else if(logical.lreflect[icv] && colvar.ss0[icv]<(cvw.lower[icv]+colvar.delta_r[icv]) && 
       colvar.ss0[icv]>cvw.lower[icv] && colvar.on[icv]) { inv_ss0[icv] = 2.*cvw.lower[icv]-colvar.ss0[icv]; wall=1;}
  }

  if(wall) {
    hills.ww[hills.n_hills] = hills.ww[hills.n_hills-1];
    fprintf(file, "%10.3f   ", mtd_data->time);
    for(icv=0;icv<ncv;icv++) {
      hills.ss0_t[hills.n_hills][icv] = inv_ss0[icv];      // new hill center
      if(colvar.on[icv]) fprintf(file, "%10.5f   ", hills.ss0_t[hills.n_hills][icv]);
    }
    for(icv=0;icv<ncv;icv++) {
      colvar.delta_s[hills.n_hills][icv] = colvar.delta_r[icv];       // new hill width
      if(colvar.on[icv]) fprintf(file, "%10.5f   ", colvar.delta_s[hills.n_hills][icv]);
    }
    wall=0;
    fprintf(file, "%10.5f\n", hills.ww[hills.n_hills]);
    hills.n_hills++;                              // hills added
  }

// INVERSION

  if (logical.do_inversion) {
    for(i_c=0;i_c<ncv;i_c++) {
      if(colvar.on[i_c]) {
        tmps=colvar.delta_r[i_c];
        tmp=colvar.ss0[i_c];
        for (j=0;j<2;j++) { // loop on 2 limits
          if (logical.invert[i_c][j]) {
            tmpd=fabs(colvar.inv_limit[i_c][j]-colvar.ss0[i_c]);
            smooth_func=1/(1+pow(tmpd/(colvar.inv_inv[i_c]*colvar.inv_ref[i_c]*tmps),10));
            if ( tmpd <= colvar.inv_ref[i_c]*tmps ) { // just mirror ...
              colvar.ss0[i_c]=2*colvar.inv_limit[i_c][j]-tmp;

              // ADD REFLECTED HILLS 
              if(hills.n_hills+10>hills.ntothills) hills_reallocate(mtd_data);
              hills.ww[hills.n_hills] = ww;
              for(icv=0;icv<ncv;icv++) hills.ss0_t[hills.n_hills][icv] = colvar.ss0[icv];   // new hill center
              for(icv=0;icv<ncv;icv++) colvar.delta_s[hills.n_hills][icv] = delta[icv];       // new hill width
              // PRINT
              if(!file) file = fopen((mtd_data->ionode?mtd_data->hilfilen:"/dev/null"), "a");

              fprintf(file, "%10.3f   ", mtd_data->time);
              for(icv=0;icv<ncv;icv++) if(colvar.on[icv]) fprintf(file, "%14.9f   ", hills.ss0_t[hills.n_hills][icv]);
              for(icv=0;icv<ncv;icv++) if(colvar.on[icv]) fprintf(file, "%14.9f   ", colvar.delta_s[hills.n_hills][icv]);
              if(logical.welltemp){
               fprintf(file, "%14.9f   %4.3f \n", hills.ww[hills.n_hills]*colvar.wfactor/(colvar.wfactor-1.0)/mtd_data->eunit,colvar.wfactor);
              } else {
               fprintf(file, "%14.9f   %4.3f \n", hills.ww[hills.n_hills]/mtd_data->eunit,0.0);
              }

              hills.n_hills++;                              // hills added

// flush all the time standalone: not a big deal...
#ifdef STANDALONE
  fclose(file);
#endif

              
              // END ADD REFLECTED HILLS    
              
              colvar.ss0[i_c]=tmp; // back to current position
            } else if ( smooth_func > 0.01 ) { // inversion condition
              hills_force();
              tmpF1=hills.Vhills; // F in current position
              colvar.ss0[i_c]=colvar.inv_limit[i_c][j];
              hills_force();
              tmpF2=hills.Vhills; // F on the border
              colvar.ss0[i_c]=2.*colvar.inv_limit[i_c][j]-tmp;
              hills_force();
              tmpF3=hills.Vhills; // F on symmetric position
              if(hills.n_hills+10>hills.ntothills) hills_reallocate(mtd_data);
              hills.ww[hills.n_hills] = (2.*tmpF2-tmpF1-tmpF3)*smooth_func;
              if (fabs(hills.ww[hills.n_hills]) > colvar.inv_maxww[i_c]*ww) {
                hills.ww[hills.n_hills] = (colvar.inv_maxww[i_c]*ww*fabs(hills.ww[hills.n_hills]))/hills.ww[hills.n_hills];
              }
    
              // ADD INVERTED HILLS

              for(icv=0;icv<ncv;icv++) hills.ss0_t[hills.n_hills][icv] = colvar.ss0[icv];   // new hill center
              for(icv=0;icv<ncv;icv++) colvar.delta_s[hills.n_hills][icv] = delta[icv];       // new hill width
              // PRINT
              if(!file) file = fopen((mtd_data->ionode?mtd_data->hilfilen:"/dev/null"), "a");

              fprintf(file, "%10.3f   ", mtd_data->time);
              for(icv=0;icv<ncv;icv++) if(colvar.on[icv]) fprintf(file, "%14.9f   ", hills.ss0_t[hills.n_hills][icv]);
              for(icv=0;icv<ncv;icv++) if(colvar.on[icv]) fprintf(file, "%14.9f   ", colvar.delta_s[hills.n_hills][icv]);
              if(logical.welltemp){
               fprintf(file, "%14.9f   %4.3f \n", hills.ww[hills.n_hills]*colvar.wfactor/(colvar.wfactor-1.0)/mtd_data->eunit,colvar.wfactor);
              } else {
               fprintf(file, "%14.9f   %4.3f \n", hills.ww[hills.n_hills]/mtd_data->eunit,0.0);
              }

              hills.n_hills++;                              // hills added

// flush all the time standalone: not a big deal...
#ifdef STANDALONE
  fclose(file);
#endif

              // END ADD INVERTED HILLS 
              
              colvar.ss0[i_c]=tmp; // back to current position
            }
          } // check inversion active
        } // loop on 2 limits
      } // check active
    } // loop on CVs, end of inversion
  }

// END INVERSION


// FLUSH FILE
  fflush(file);
  if(logical.do_walkers){
// we close it for multiple walkers calculation, so as to allow other processes to read it
    fclose(file);
    file=NULL;
  }

}

// This routine add the gaussian hills
void PREFIX hills_add(struct mtd_data_s *mtd_data)

{
  int icv,irep;
  real* all_ww;
  real** all_ss;
  real** all_delta;
  real  this_ww;
  real  this_ss[nconst_max];
  real  this_delta[nconst_max];
  int ineighbour,distance;
  int nrep,ncv;
  int myrep;
  static int last_hill_at_this_step;
  static int first=1;

  if(first) last_hill_at_this_step=colvar.it;
  first=0;

  ncv=colvar.nconst;

  if(hills.max_height>0.0) hills.wwr=hills.rate*(colvar.it-last_hill_at_this_step)*mtd_data->dt;

  // ADD THE HILL
  // new hill height 
  this_ww = hills.wwr;
  if(logical.welltemp) {
  /* since hills are added after the calculation of the bias potential, we can reuse
     the stored Vhills to decide the hills height*/
    this_ww *= exp(-hills.Vhills/(mtd_data->boltz*(colvar.wfactor-1.0)*colvar.simtemp));
  }
  // JFD>
  if(logical.transition_tempering) {
    this_ww *= exp(-transition_bias_ND()/(mtd_data->boltz*(colvar.ttfactor-1.0)*colvar.simtemp));
  }
  // <JFD
  // ADW>
  if(logical.target_distribution) {
    this_ww /= exp(mtd_data->boltz *colvar.simtemp * grid_getstuff(&target_grid, colvar.ss0,  NULL));
    if(logical.welltemp)//to prevent very large hills
      this_ww = fmin(hills.wwr, this_ww);
#ifdef OUTPUT_HILL_MULTIPLIERS    
    printf("Hill = fmin(%f, %f * %f * %f =  %f)\n",this_ww, this_ww, exp(-hills.Vhills/(mtd_data->boltz*(colvar.wfactor-1.0)*colvar.simtemp)), 1. / exp(mtd_data->boltz * colvar.simtemp * grid_getstuff(&target_grid, colvar.ss0,  NULL)), hills.wwr * exp(-hills.Vhills/(mtd_data->boltz*(colvar.wfactor-1.0)*colvar.simtemp)) / exp(mtd_data->boltz * colvar.simtemp * grid_getstuff(&target_grid, colvar.ss0,  NULL)));
#endif
  }
  if(logical.mcgdp_hills) {
    //check if the colvar is within the boundaries
    for (icv=0; icv<ncv; icv++) if(colvar.on[icv]) {
	if (hills.mcgdp_reshape_flag[icv] == 1) {
	  if(colvar.ss0[icv] < hills.hill_lower_bounds[icv] ||
	     colvar.ss0[icv] > hills.hill_upper_bounds[icv]) {
	    //we aren't, make hill_height 0
	    this_ww = 0;
	    break;
	  }
	    
	}
      }
  }
  // <ADW

  
  for(icv=0;icv<ncv;icv++) this_ss[icv]    = colvar.ss0[icv];    	// new hill center
  for(icv=0;icv<ncv;icv++) this_delta[icv] = colvar.delta_r[icv];       // new hill width

  if(hills.max_height>0.0 && this_ww<hills.max_height && (colvar.it-last_hill_at_this_step)<hills.max_stride) return;
  last_hill_at_this_step=colvar.it;

  
// add hill to the array
  if(! (logical.remd && colvar.ptmetad_neighbours>0) ) hills_push(mtd_data,this_ww,this_ss,this_delta);

// NEIGHBOURS HILLS for parallel tempering
  if(logical.remd && colvar.ptmetad_neighbours>0){
    nrep=mtd_data->nrepl;
    myrep=mtd_data->repl;
    all_ww=float_1d_array_alloc(nrep);
    all_ss=float_2d_array_alloc(nrep,ncv);
    all_delta=float_2d_array_alloc(nrep,ncv);
    for(irep=0;irep<nrep;irep++)                          all_ww[irep]=0.0;
    for(irep=0;irep<nrep;irep++) for(icv=0;icv<ncv;icv++) all_ss[irep][icv]=0.0;
    for(irep=0;irep<nrep;irep++) for(icv=0;icv<ncv;icv++) all_delta[irep][icv]=0.0;
// first broadcast on the master node of each replica
    if(mtd_data->ionode){
      all_ww[myrep]=this_ww;
      for(icv=0;icv<ncv;icv++) all_ss[myrep][icv]=this_ss[icv];
      for(icv=0;icv<ncv;icv++) all_delta[myrep][icv]=this_delta[icv];
      plumed_intersum(mtd_data,nrep, & all_ww[0]);
      plumed_intersum(mtd_data,nrep*ncv, & all_ss[0][0]);
      plumed_intersum(mtd_data,nrep*ncv, & all_delta[0][0]);
    }
// then broadcast inside each replica
    plumed_sum(mtd_data,nrep, & all_ww[0]);
    plumed_sum(mtd_data,nrep*ncv, & all_ss[0][0]);
    plumed_sum(mtd_data,nrep*ncv, & all_delta[0][0]);
    for(ineighbour=0;ineighbour<nrep;ineighbour++){
      distance=myrep-ineighbour;
      if(distance<0) distance=-distance;
      if(distance<=colvar.ptmetad_neighbours){
        this_ww=all_ww[ineighbour]*exp(-0.5*distance*distance/(colvar.ptmetad_sigma*colvar.ptmetad_sigma));;
        for(icv=0;icv<ncv;icv++)this_ss[icv]=all_ss[ineighbour][icv];
        for(icv=0;icv<ncv;icv++)this_delta[icv]=all_delta[ineighbour][icv];
        hills_push(mtd_data,this_ww,this_ss,this_delta);
      }
    }
    free_1dr_array_alloc(all_ww);
    free_2dr_array_alloc(all_ss,nrep);
    free_2dr_array_alloc(all_delta,nrep);
  }

// After hill addition, update mtd_data.Vbias
  hills.Vhills=hills_engine(colvar.ss0,colvar.ff_hills);
}

//-----------------------------------------------------------------------------------------------

real PREFIX hills_engine_dp(int ih,real* ss0,real* dp){
  int icv,ncv;
  real dp2,diff;
  ncv=colvar.nconst;
  dp2 = 0.;
  for(icv=0;icv<ncv;icv++)if(colvar.on[icv]){
    dp[icv] = (ss0[icv]-hills.ss0_t[ih][icv])/colvar.delta_s[ih][icv];       // (cv(now)-hil(ih))/sigma
    if(colvar.type_s[icv]==5 || (  colvar.type_s[icv]==34 &&  colvar.type[icv]==2 )) {
      diff = ss0[icv]-hills.ss0_t[ih][icv];
      if(diff>M_PI) diff-=2.*M_PI;
      else if(diff<-M_PI) diff+=2.*M_PI;
      dp[icv] = diff/colvar.delta_s[ih][icv];
    }
    dp2 += dp[icv]*dp[icv];                  // sum dp*dp
  }
  dp2 = 0.5*dp2;
  return dp2;
};

real PREFIX hills_engine(real* ss0,real* force){
/* this routine calculate hills forces and energy in an arbitrary point.
   let us try to use it anywhere it is needed, so as to avoid errors and double coding.
   in particular, ADD HERE VARIABLES NEEDING PBC
   In GROMACS, DL_POLY, and AMBER it is parallel, and should be called
   with the same arguments by all the processes belongin to a replica
*/

  int dp2index;
  real dp2, dp[nconst_max], VhillsLast;
  real Vbias;
  // JFD>
  // Temporaries for calculating McGovern-de Pablo boundary-consistent
  // hills and forces.
  real mcgdp_VHillDenom;
  real lbound_exp_argument, ubound_exp_argument;
  real lbound_gaussian, ubound_gaussian;
  real mcdgp_force_correction;
  int erf_index, lbound_exp_index, ubound_exp_index;
  // <JFD

  int nh,ncv;
  int ih,icv;

  int npe,rank;

  int nhstart; /* first hill belonging to this process */
  int nhstride; /* stride for hills belonging to this process */

  nh=hills.n_hills;
  ncv=colvar.nconst;

  if(logical.parallel_hills){
    npe=plumed_comm_size(&mtd_data);
    rank=plumed_comm_rank(&mtd_data);
  }else{
    npe=1;
    rank=0;
  };

  if(logical.do_grid){
/*
   when grid are used, it is better NOT to parallelize on hill index.
   in fact, since grid are based on a non-linear spline, splitting the hills among different processors
   leads to slighlty different forces, which means non reproducibility in parallel calculations.
   moreover, the advantage of this parallelization is only when restarting (usually one hill at a time is added)
*/
    nhstart=0;
    nhstride=1;
  }else{
    nhstart=rank;
    nhstride=npe;
  };

  Vbias=0.0;
  if(force) for(icv=0;icv<ncv;icv++) force[icv]=0.0;

  //ADW: The loop below is way too slow as the number of hills is
  //increased on a grid. I'm separating it for the different cases.
  if(logical.do_grid && !logical.debug_grid) { 
    for(ih=bias_grid.nhills; ih < nh && ih < hills.read; ih+=nhstride){
      grid_addhills(&bias_grid,hills.ww[ih],hills.ss0_t[ih],colvar.delta_s[ih],rank,npe); 
    }
    nhstart = hills.read;
  }
  
// This loop is parallelized
  // if logical.debug_grid is set, the actual force is calculated with the hills, but
  // in a debug file we write the actual and the grid force and energy

  for(ih=nhstart;ih<nh;ih+=nhstride){
    dp2 = hills_engine_dp(ih, ss0, dp);
    if(dp2 < DP2CUTOFF){
      dp2index =  dp2 * GTAB / DP2CUTOFF;
      VhillsLast = hills.ww[ih] * hills.exp[dp2index];
      // JFD>
      // The McGovern-de Pablo boundary consistent hills require dividing
      // the usual hill Gaussian by an error function term representing
      // the integral of that Gaussian over the allowed interval.
      if (logical.mcgdp_hills) {
	mcgdp_VHillDenom = 1.0;
	for (icv=0; icv<ncv; icv++) if(colvar.on[icv]) {
            if (hills.mcgdp_reshape_flag[icv] == 1) {
              erf_index = (ss0[icv] - hills.hill_lower_bounds[icv])/(hills.hill_upper_bounds[icv]-hills.hill_lower_bounds[icv]) * GTAB;
              mcgdp_VHillDenom *= hills.erf[erf_index][icv] / 2.0;
            }
          }
	VhillsLast /= mcgdp_VHillDenom;
      }
      // <JFD
      Vbias += VhillsLast;
      if(force) for(icv=0;icv<ncv;icv++) if(colvar.on[icv]) {
	    if(logical.interval[icv]) {
	      if((ss0[icv]> cvint.lower_limit[icv] && ss0[icv]<cvint.upper_limit[icv])) {
		force[icv] += dp[icv] / colvar.delta_s[ih][icv] * VhillsLast;  // -dU/dCV
	      }
	      // JFD>
	      // Applying the product rule to McGovern-de Pablo hill, taking
	      // derivatives of the numerator and denominator into separate
	      // terms, results in one normal-looking but rescaled Gaussian 
	      // force and one more complex boundary correction term.
	    } else if (logical.mcgdp_hills && hills.mcgdp_reshape_flag[icv]) {
	      // Compute how large the hill is where it hits the lower bound
	      lbound_exp_argument = (ss0[icv] - hills.hill_lower_bounds[icv]) * (ss0[icv] - hills.hill_lower_bounds[icv])/(2 * colvar.delta_s[ih][icv] * colvar.delta_s[ih][icv]);
	      if (lbound_exp_argument < DP2CUTOFF) {
		lbound_exp_index = lbound_exp_argument * GTAB / DP2CUTOFF;
		lbound_gaussian = hills.exp[lbound_exp_index];
	      } else {
		lbound_gaussian = 0;
	      }
	      // Compute how large the hill is where it hits the upper bound
	      ubound_exp_argument = (hills.hill_upper_bounds[icv] - ss0[icv]) * (hills.hill_upper_bounds[icv] - ss0[icv])/(2 * colvar.delta_s[ih][icv] * colvar.delta_s[ih][icv]);
	      if (ubound_exp_argument < DP2CUTOFF) {
		ubound_exp_index = ubound_exp_argument * GTAB / DP2CUTOFF;
		ubound_gaussian = hills.exp[ubound_exp_index];
	      } else {
		ubound_gaussian = 0;
	      }
	      // Compute the term corresponding to taking the derivative of the
            // hill function's denominator
	      if (lbound_gaussian > 0 || ubound_gaussian > 0) {
		mcdgp_force_correction = M_sqrt2oPI * (lbound_gaussian - ubound_gaussian) * VhillsLast / hills.erf[erf_index][icv] * colvar.delta_s[ih][icv];
	      } else {
		mcdgp_force_correction = 0;
	      }
	      force[icv] += dp[icv] / colvar.delta_s[ih][icv] * VhillsLast + mcdgp_force_correction;
	      // <JFD
	    } else {
	      force[icv] += dp[icv] / colvar.delta_s[ih][icv] * VhillsLast;  // -dU/dCV
	    }
	  }
    }
    
  }


  if(logical.do_grid) {
    bias_grid.nhills = hills.read;
    if(!logical.debug_grid){
      Vbias+=grid_getstuff(&bias_grid,ss0,force);
    } else {
      static FILE *file=NULL;
      if(!file) {
        char buf[1024];
        sprintf(buf,"DEBUG_GRID%i+%i",mtd_data.repl,plumed_comm_rank(&mtd_data));
        fprintf(stderr,"%s\n",buf);
        file=fopen(buf,"w");
        
      }
      real VbiasG;
      real forceG [nconst_max];
      int i;
      for(i=0;i<ncv;i++)  forceG[i]=0;
      VbiasG=grid_getstuff(&bias_grid,ss0,forceG);
      for(i=0;i<ncv;i++) fprintf(file," %f",ss0[i]);
      fprintf(file,"   %f %f    ",Vbias,VbiasG);
      if(force) for(i=0;i<ncv;i++) fprintf(file," %f %f",force[i],forceG[i]);
      fprintf(file,"\n");
    }
  }

/* when using grid, each process knows the whole grid, so there is no need to reduce */
  if(logical.parallel_hills && ! logical.do_grid){
    if(force) plumed_sum(&mtd_data,ncv, & force[0]);
    plumed_sum(&mtd_data,1, & Vbias);
  }
  
  return Vbias;
};

// This routine calculate the hills contribution
void PREFIX hills_force()
{
  hills.Vhills=hills_engine(colvar.ss0,colvar.ff_hills);
}

//-----------------------------------------------------------------------------------------------

// This routine read the HILLS file
void PREFIX read_hills(struct  mtd_data_s *mtd_data, int restart, int first_read)
{
  double dummy;
  real  old_factor;
  int i, j, nactive, iw;
  long int line;
  FILE *file;
  char *str, stringa[64000];
  char walkfilen[64000];

  if(restart) hills.first_read=0;
  
  /****************************** UPDATING *************************/

  for(iw=0;iw<hills.nwalkers;iw++){                                                             // cycle on walkers

   file=NULL;
  
   if(logical.do_walkers) {
    sprintf(walkfilen, "%s.%i", mtd_data->basefilen, iw);  
   } else {
    sprintf(walkfilen, "%s", mtd_data->hilfilen); 
   }

   file = fopen(walkfilen, "r");
   if(!file && logical.do_walkers) continue;
   if(!file && !logical.do_walkers) {										// not exist
           char buf[1024];
           sprintf(buf,"Cannot read HILLS file :: %s",mtd_data->hilfilen);
           plumed_error(buf);
   }
   fflush(file); 

   if(!restart && !first_read)fsetpos(file,&(hills.line_counter[iw]));
 
   line = hills.read;
 
   int active_to_ind[nconst_max];
   nactive=0.; //number of active cvs 
   for(i=0;i<colvar.nconst;i++){
      if(colvar.on[i]){active_to_ind[nactive]=i;nactive++;}
   }
 
   while(1){											// read cycle
     str = fgets(stringa, 64000, file);
     if(str == NULL) break;
 
     // reallocate if needed during reading
     if(line+10>hills.ntothills) hills_reallocate(mtd_data);  
 
     j=0;
     str =(char *) strtok(stringa," \t");
     // skip header line
     if(strcmp(str,"#!")==0||strcmp(str,"#")==0){ continue; }
     while (str != NULL)
     {
       if( j>0 && j <=nactive  ) { 
           i=active_to_ind[j-1];
           sscanf(str, "%lf", &dummy);
           hills.ss0_t[line][i] = (real) dummy; 
       //    printf("POS %d   %f ",i,  hills.ss0_t[line][i] );
       }
       else if( j>nactive && j<= 2*nactive ) { 
           i=active_to_ind[j-nactive-1];
           sscanf(str, "%lf", &dummy);
           colvar.delta_s[line][i] = (real) dummy;							// read the hills dimension
       //    printf("DELTA %d   %f ",i, colvar.delta_s[line][i]);
       }
       else if( j==2*nactive+1 ) { 
           sscanf(str, "%lf", &dummy);
           hills.ww[line] = (real) dummy * mtd_data->eunit;                                                              // read the hills height
       //   printf("WW   %f \n", hills.ww[line]);
       }
       else if( j==2*nactive+2 ) {
           sscanf(str, "%lf", &dummy);
           old_factor = (real) dummy;
           if(logical.welltemp && !logical.read_old_bf) hills.ww[line] = hills.ww[line] * (colvar.wfactor-1.0) / colvar.wfactor;
           if(logical.welltemp && logical.read_old_bf){
            if(old_factor<1.0) plumed_error("Restarting from a non-welltempered metadynamics. Please remove READ_OLD_BF."); 
            hills.ww[line] = hills.ww[line] * (old_factor-1.0) / old_factor;
//            printf("BF   %f \n", old_factor);
           }
       }
       str =(char *) strtok(NULL, " \t");
       j++;
     }
     line++;
   }
   
   fgetpos(file,&(hills.line_counter[iw]));
 
   fclose(file);
   hills.n_hills = line;
 
   if(restart){
     fprintf(mtd_data->fplog, "|- RESTARTING HILLS: TOT %li HILLS read from %s\n",hills.n_hills-hills.read,walkfilen);
   }else{
     if(hills.n_hills-hills.read>0) 
      fprintf(mtd_data->fplog, "|- UPDATING HILLS: from %li to %li TOT %li HILLS read from %s \n", hills.read,hills.n_hills-1,hills.n_hills-hills.read,walkfilen);
   }
   hills.read=hills.n_hills;

  } // end cycle on walkers
  
//  for(i=hills.read;i<hills.n_hills;i++) fprintf(mtd_data->fplog, "UPDATING # %i HILLS %f %f \n",i,hills.ss0_t[i][0],hills.ss0_t[i][1]);
//  for(i=0;i<hills.n_hills;i++) fprintf(mtd_data->fplog, "AFTER UPDATE # %i HILLS %f %f \n",i,hills.ss0_t[i][0],hills.ss0_t[i][1]);
  fprintf(mtd_data->fplog,"\n");
  fflush(mtd_data->fplog);
}

//------------------------------------------------------------------------------------------

void PREFIX hills_reallocate( struct mtd_data_s *mtd_data) { 
      int i,j;
      long int oldtot;
      real  *ww_tmp;
      real  **ss0_t_tmp;
      real  **delta_s_tmp; 

      oldtot=hills.ntothills;
      hills.ntothills+=STACKDIM ;

      fprintf(mtd_data->fplog,"HILLS REALLOCATION--OLD DIMENSION: %li\n",oldtot);
      // save pointers to old arrays
      ww_tmp        = hills.ww;
      ss0_t_tmp     = hills.ss0_t;
      delta_s_tmp   = colvar.delta_s;

      // reallocate to the new stackdim
      hills.ww        = float_1d_array_alloc(hills.ntothills);
      hills.ss0_t     = float_2d_array_alloc(hills.ntothills,colvar.nconst); 
      colvar.delta_s  = float_2d_array_alloc(hills.ntothills,colvar.nconst);   

      // copy old arrays to new ones
      // WE COPY THE FULL ARRAY, NOT JUST hills.n_hills.
      // this is safer, since in reallocation during reading n_hills only set at the end
      for (i=0;i<oldtot;i++){
            hills.ww[i]=ww_tmp[i];  
            for (j=0;j<colvar.nconst;j++) hills.ss0_t[i][j]=ss0_t_tmp[i][j]; 
            for (j=0;j<colvar.nconst;j++) colvar.delta_s[i][j]=delta_s_tmp[i][j]; 
      } 

      // free old arrays
      free_1dr_array_alloc(ww_tmp);
      free_2dr_array_alloc(ss0_t_tmp,oldtot);
      free_2dr_array_alloc(delta_s_tmp,oldtot);

      fprintf(mtd_data->fplog,"HILLS REALLOCATION--NEW DIMENSION: %li\n",hills.ntothills);
}
//-------------------------------------------------------------------------------------------
// Calculate total grid dimension, allocate and initialize
void PREFIX grid_initialize(struct grid_s *grid)
{

 int i, j, ncv;
 
 grid->one2multi      = NULL;
 grid->one2multi_full = NULL;
 grid->size           = 1;
 ncv                  = grid->ncv;
 
 for(i=0;i<ncv;i++) {
  grid->lbox[i]  = grid->max[i] - grid->min[i];  
  grid->dx[i]    = grid->lbox[i] / grid->bin[i];
  if(grid->period[i]==0){
   grid->max[i]  += grid->dx[i];
   grid->lbox[i] += grid->dx[i];
   grid->bin[i]  += 1;
  }
  grid->size    *= grid->bin[i];
 }

// allocation
 grid->pot   = float_1d_array_alloc(grid->size);  
 grid->force = float_2d_array_alloc(grid->size,ncv);
 grid->one2multi_full=int_2d_array_alloc(grid->size,grid->ncv); 
// filling one2multi_full
 grid_create_one2multi(grid->one2multi_full, grid->size, grid->ncv, grid->bin); 

 grid->mem = (ncv+1)*grid->size*sizeof(real)/pow(1024.0,2);
//fprintf(stderr,"FFF %i\n",plumed_comm_rank(&mtd_data));
 fprintf(mtd_data.fplog,"|- GRID MEMORY USAGE :: %6.2f MB \n",grid->mem);

 for(j=0;j<grid->size;j++) {
  grid->pot[j] = 0. ;
  for(i=0;i<ncv;i++) grid->force[j][i] = 0. ;
 }

}
//-------------------------------------------------------------------------------------------
// Calculate reduced grid dimension. This routine is called once unless
// the Gaussian sigma (colvar.delta) is modified during the simulation.
void PREFIX grid_resize_minigrid(struct grid_s *grid, real* delta, real cutoff)
{

 real mem;
 int  i, ncv;

// store cutoff for later use
 grid->cutoff   = cutoff;
 ncv            = grid->ncv;
 grid->minisize = 1;

 for(i=0;i<ncv;i++) {
  grid->minilbox[i] = sqrt(2.*cutoff)*delta[grid->index[i]];      // this is HALF the side of minibox
  grid->minibin[i]  = floor(2.*grid->minilbox[i]/grid->dx[i])+1;
  grid->minisize   *= grid->minibin[i];
 } 

 fprintf(mtd_data.fplog,"|- UPDATING REDUCED GRID: SIZE %d pts on %d ",grid->minisize,grid->size);
 fprintf(mtd_data.fplog," DIMENSION "); for(i=0;i<ncv-1;i++) fprintf(mtd_data.fplog," %d x",grid->minibin[i]); 
 fprintf(mtd_data.fplog," %d \n",grid->minibin[ncv-1]);

// deallocate if previously allocated
 if(grid->one2multi) free_2di_array_alloc(grid->one2multi,grid->minisize);
// new allocation
 grid->one2multi = int_2d_array_alloc(grid->minisize,grid->ncv);
 mem = grid->minisize*grid->ncv*sizeof(int)/pow(1024.0,2);
 fprintf(mtd_data.fplog,"|- GRID MEMORY USAGE :: %6.2f MB \n", grid->mem+mem);

// create one2multi vector
 grid_create_one2multi(grid->one2multi, grid->minisize, grid->ncv, grid->minibin);

}
//-------------------------------------------------------------------------------------------
// Allocate and create one2multi/one2multi_full vector. 
// One2multi is needed to work efficiently on the reduced grid:
// the routine is called once for the minigrid unless the HILLS delta is modified during the simulation.
// One2multi_full is used when writing or restarting a simulation from a GRID on file
void PREFIX grid_create_one2multi(int **one2multi, int size, int ncv, int *bin)
{

 int i, j, k, tmpgrid, index1d;

 if(ncv==1) {
  for(i=0;i<size;i++) one2multi[i][0] = i;
 } else {
  for(i=0;i<size;i++){
   index1d=i;
   for(j=ncv-1;j>=0;j--){
    tmpgrid = 1;
    for(k=0;k<j;k++) tmpgrid*=bin[k];
    one2multi[i][j] = index1d/tmpgrid;
    index1d = index1d%tmpgrid;
    if(index1d==-1) {
     one2multi[i][j] -= 1;
     index1d=tmpgrid;
    }
   }
  } 
 }

}
//-------------------------------------------------------------------------------------------
// add a hills on the grid. Update potential and force
void PREFIX grid_addhills(struct grid_s *grid, real ww, real* ss, real* delta,int rank,int npe)
{

  int   i, j, ncv, flag;
  real *xx, *dp, dp2, expo;
  int  *index_nd, index_1d, dp2index;
  int *index_1d_para;
  real *pot_for_para;

  // JFD>
  // Temporaries for the Dicksonian tempering rule
  real bias_at_hill_center;
  real bias_at_grid_point;
  real *force_at_hill_center;
  real *force_at_grid_point;
  real dicksonian_biasing_factor;
  real inv_tempering_temp = 1.0 / (mtd_data.boltz*(colvar.wfactor-1.0)*colvar.simtemp);
  // Temporaries for calculating McGovern-de Pablo boundary-consistent
  // hills and forces.
  real mcgdp_VHillDenom;
  real lbound_exp_argument, ubound_exp_argument;
  real lbound_gaussian, ubound_gaussian;
  real mcdgp_force_correction;
  int erf_index, lbound_exp_index, ubound_exp_index;
  // <JFD

  ncv  = grid->ncv;

  // allocate temp array
  xx = float_1d_array_alloc(ncv);
  dp = float_1d_array_alloc(ncv);
  index_nd = int_1d_array_alloc(ncv); 
  force_at_hill_center = float_1d_array_alloc(ncv);
  force_at_grid_point = float_1d_array_alloc(ncv);

// preliminary checks
// 1) if the HILLS center is inside the grid
// 2) if the GRID bin size is too large
// 3) if delta is changed from previous call
  flag = 0.;
  for(j = 0; j < ncv; j++) {
    if((ss[grid->index[j]] < grid->min[j] || ss[grid->index[j]] >= grid->max[j]) && !grid->period[j]) {
      plumed_error("HILLS outside GRID. Please increase GRID size."); 
    }
    if(grid->dx[j] > delta[grid->index[j]] / 2.0) plumed_error("GRID bin size is too large compared to HILLS sigma."); 
    if(fabs((grid->oldelta[j] - delta[grid->index[j]]) / delta[grid->index[j]]) > 0.05) flag = 1; 
  }
  // recalculate the dimension of the reduced grid if delta is changed
  if(flag == 1) {
    grid_resize_minigrid(grid, delta, DP2CUTOFF);
    for(j = 0; j < ncv; j++) grid->oldelta[j] = delta[grid->index[j]];
  }

  // temporary array for parallel computation
  index_1d_para = int_1d_array_alloc(grid->minisize);
  pot_for_para = float_1d_array_alloc(grid->minisize * (1 + ncv));
  for(i = 0; i < grid->minisize; i++) index_1d_para[i] = 0;
  for(i = 0; i < grid->minisize * (1 + ncv); i++) pot_for_para[i] = 0.0;

  // add HILL to the points belonging to the reduced GRID

  // JFD>  
  if(logical.dicksonian_tempering) {
    bias_at_hill_center = grid_getstuff(grid, ss, force_at_hill_center);
  }
  // <JFD

  for(i = rank; i < grid->minisize; i += npe) {

    index_1d_para[i] = -1; // it means "no force on this point"

    flag = 0;
    for(j = 0; j < ncv; j++) {
      xx[j] = ss[grid->index[j]] - grid->minilbox[j] + grid->dx[j] * grid->one2multi[i][j];
      if(grid->period[j]) xx[j] -= grid->lbox[j] * rint(xx[j]/grid->lbox[j]);
      index_nd[j] = floor((xx[j]-grid->min[j])/grid->dx[j]);
      //with single precision, it is possible that xx - min = 2 * min if xx[j] is very close to min
      if(index_nd[j]<0 || index_nd[j] >=grid->bin[j])  flag=1;
    }
    if(flag == 1) continue; // out of grid 

    // from multidimensional index to mono
    index_1d = grid_multi2one(grid, index_nd);

    // add the gaussian on the GRID
    dp2 = 0.;
    for(j = 0; j < ncv; j++) {
      xx[j] = grid->min[j] + grid->dx[j] * index_nd[j];
      dp[j] = xx[j] - ss[grid->index[j]];
      if(grid->period[j]) dp[j] -= grid->lbox[j] * rint(dp[j]/grid->lbox[j]); 
      dp[j] /= delta[grid->index[j]];
      dp2 += dp[j] * dp[j]; 
    }
    dp2 *= 0.5;

    if(dp2 < grid->cutoff){    
      dp2index =  dp2 * GTAB / DP2CUTOFF;
      expo     = ww * hills.exp[dp2index];
      // JFD>
      // Dicksonian tempering is tempering using the bias at each point to be biased rather than
      // the hill center. That can be done within this function by computing the ratio of the
      // usual well-tempering (which has already been applied) to the Dicksonian well-tempering
      // and then applying that ratio to the hill value that is calculated normally.
      if(logical.dicksonian_tempering) {
        bias_at_grid_point = grid_getstuff(grid, xx, force_at_grid_point);
        dicksonian_biasing_factor = exp(-inv_tempering_temp * (bias_at_grid_point - bias_at_hill_center));
        expo *= dicksonian_biasing_factor;
      }
      // The McGovern-de Pablo boundary consistent hill requires dividing
      // the usual hill Gaussian by an error function term representing
      // the integral of that Gaussian over the allowed interval.
      if (logical.mcgdp_hills) {
        mcgdp_VHillDenom = 1.0;
        for (j = 0; j < ncv; j++) if(colvar.on[j]) {
          if (hills.mcgdp_reshape_flag[j] == 1) {
	    if(xx[j] > hills.hill_lower_bounds[j] && xx[j] < hills.hill_upper_bounds[j]) {
	      erf_index = (xx[j] - hills.hill_lower_bounds[j])/(hills.hill_upper_bounds[j]-hills.hill_lower_bounds[j]) * GTAB;
	      mcgdp_VHillDenom *= hills.erf[erf_index][j] / 2.0;
	    }
          }
        }
        expo /= mcgdp_VHillDenom;
      }
      // <JFD
      
      // Add grid bias potential value
      pot_for_para[i * (ncv + 1)] = expo;
      // Add grid bias force vector
      for(j = 0; j < ncv; j++) pot_for_para[i * (ncv + 1) + 1 + j] = dp[j] / delta[grid->index[j]] * expo;
      
      // JFD>
      // Because Dicksonian tempering is location-sensitive, the forces to add depend on the current
      // bias force at the point.
      if(logical.dicksonian_tempering) {
        for(j = 0; j < ncv; j++) {
          pot_for_para[i * (ncv + 1) + 1 + j] += inv_tempering_temp * force_at_grid_point[j] * expo;
        }
      }
      // Applying the product rule to McGovern-de Pablo hill, taking
      // derivatives of the numerator and denominator into separate
      // terms, results in one normal-looking but rescaled Gaussian 
      // force and one more complex boundary correction term. The 
      // numerator term has already been added to the force at this point,
      // so this only calculates the denominator term.
      if (logical.mcgdp_hills) {
        for(j = 0; j < ncv; j++) {
          if (hills.mcgdp_reshape_flag[j] == 1) {
	    //check if we're actually within the McGovern-de Pablo bounds
	    if(xx[j] <= hills.hill_lower_bounds[j] || xx[j] >= hills.hill_upper_bounds[j])
	      continue;
	    
            // Compute how large the hill is where it hits the lower bound
            lbound_exp_argument = (xx[j] - hills.hill_lower_bounds[j]) * (xx[j] - hills.hill_lower_bounds[j])/(2 * delta[grid->index[j]] * delta[grid->index[j]]);
            if (lbound_exp_argument < DP2CUTOFF) {
              lbound_exp_index = lbound_exp_argument * GTAB / DP2CUTOFF;
              lbound_gaussian = hills.exp[lbound_exp_index];
            } else {
              lbound_gaussian = 0;
            }
            // Compute how large the hill is where it hits the upper bound
            ubound_exp_argument = (hills.hill_upper_bounds[j] - xx[j]) * (hills.hill_upper_bounds[j] - xx[j])/(2 * delta[grid->index[j]] * delta[grid->index[j]]);
            if (ubound_exp_argument < DP2CUTOFF) {
              ubound_exp_index = ubound_exp_argument * GTAB / DP2CUTOFF;
              ubound_gaussian = hills.exp[ubound_exp_index];
            } else {
              ubound_gaussian = 0;
            }
            // Compute the term corresponding to taking the derivative of the
            // hill function's denominator
            if (lbound_gaussian > 0 || ubound_gaussian > 0) {
              mcdgp_force_correction = M_sqrt2oPI * (lbound_gaussian - ubound_gaussian) * expo / hills.erf[erf_index][j] * delta[grid->index[j]];
            } else {
              mcdgp_force_correction = 0;
            }
            pot_for_para[i * (ncv + 1) + 1 + j] += mcdgp_force_correction;
          }
        }
      }
      // <JFD
      index_1d_para[i]=index_1d;
    }
  }

  if(npe>1){ 
    plumed_sum (&mtd_data,grid->minisize*(ncv+1),pot_for_para);
    plumed_sumi(&mtd_data,grid->minisize,index_1d_para);
  }

  for(i=0;i<grid->minisize;i++) {
    if(index_1d_para[i]<0) continue;
    grid->pot[index_1d_para[i]]+=pot_for_para[i*(ncv+1)];
    for(j=0;j<ncv;j++) grid->force[index_1d_para[i]][j] += pot_for_para[i*(ncv+1)+1+j];
  }

  // deallocation
  free_1dr_array_alloc(dp);
  free_1dr_array_alloc(xx);
  free_1di_array_alloc(index_nd);
  free_1dr_array_alloc(pot_for_para);
  free_1di_array_alloc(index_1d_para);
}

//-------------------------------------------------------------------------------------------
// from multidimensional index to mono dimensional
int PREFIX grid_multi2one(struct grid_s *grid, int* index_nd)
{
 int i, j, index, tmpgrid;

 index = index_nd[0];
 
 for(i=1;i<grid->ncv;i++) {
  tmpgrid = 1;
  for(j=0;j<i;j++) tmpgrid *= grid->bin[j];
  index += index_nd[i]*tmpgrid;
 }
  
 return index;

}
//-------------------------------------------------------------------------------------------
// this routine returns the bias in a certain point and optionally the forces
real PREFIX grid_getstuff(struct grid_s *grid, real* ss0, real* force)
{

 real *xx;
 int  j, *index_nd, index_1d, ncv;
 real Vbias;

 ncv  = grid->ncv;

// allocate temp array
 xx       = float_1d_array_alloc(ncv);
 index_nd = int_1d_array_alloc(ncv);

// first check if the point is inside the GRID
 for(j=0;j<ncv;j++) {
  xx[j] = ss0[grid->index[j]];
  if((xx[j]<grid->min[j] || xx[j]>=grid->max[j]) && !grid->period[j]) {
      fprintf(stderr, "Bad grid value = %f\n", xx[grid->index[j]]);	
      plumed_error("You are outside the GRID!. Please increase GRID size.");
  }
  index_nd[j] = floor((xx[j]-grid->min[j])/grid->dx[j]);
  if(grid->period[j]) index_nd[j] -= grid->bin[j] * int_floor((real) index_nd[j]/grid->bin[j]);
 }

// from multidimensional index to mono
 index_1d=grid_multi2one(grid,index_nd);

 if(!logical.donot_spline){
   real f;
   real where[nconst_max];
   int  stride[nconst_max];
   real der[nconst_max];
   for(j=0;j<ncv;j++) where[j]=xx[j]-grid->min[j]-index_nd[j]*grid->dx[j];
   stride[0]=1;
   for(j=1;j<ncv;j++) stride[j]=stride[j-1]*grid->bin[j-1];
   for(j=0;j<ncv;j++) if(grid->period[j]  && index_nd[j]==grid->bin[j]-1) stride[j]*=(1-grid->bin[j]);
   for(j=0;j<ncv;j++) if(!grid->period[j] && index_nd[j]==grid->bin[j]-1) plumed_error("You are outside the GRID!. Please increase GRID size.");
   f=spline(ncv,grid->dx,where,& grid->pot[index_1d],&grid->force[index_1d][0],stride,der);
   Vbias=f;
   if(force) for(j=0;j<ncv;j++) force[grid->index[j]] -=der[j];
 } else {
// getting BIAS and FORCE
   Vbias = grid->pot[index_1d];
   if(force) for(j=0;j<ncv;j++) force[grid->index[j]] += grid->force[index_1d][j];
 }

// free
 free_1dr_array_alloc(xx);
 free_1di_array_alloc(index_nd);

 return Vbias;
}
//-------------------------------------------------------------------------------------------
// write GRID to file
void  PREFIX grid_write_tofile(struct grid_s *grid)
{

 int      i, j, *index_nd, *bin;
 real   *xx, *max;
 FILE     *file=NULL;

// Open grid file for writing 
 file = fopen(grid->w_file, "w");

// Allocate stuff
 xx       = float_1d_array_alloc(grid->ncv);
 max      = float_1d_array_alloc(grid->ncv);
 index_nd = int_1d_array_alloc(grid->ncv);
 bin      = int_1d_array_alloc(grid->ncv); 

// if the cv is not periodic we have to subtract 
// in output 1 to bin size and dx to max
// to respect our convention
 for(i=0;i<grid->ncv;i++){
  if(grid->period[i]==0){
   bin[i]=grid->bin[i]-1;
   max[i]=grid->max[i]-grid->dx[i];
  } else {
   bin[i]=grid->bin[i];
   max[i]=grid->max[i];
  }
 }

// HEADER
 fprintf(file,"#! FORCE 1\n");
 fprintf(file,"#! NVAR %d\n",grid->ncv);
 fprintf(file,"#! TYPE"); for(i=0;i<grid->ncv;i++) fprintf(file," %d",  colvar.type_s[grid->index[i]]); fprintf(file,"\n");
 fprintf(file,"#! BIN");  for(i=0;i<grid->ncv;i++) fprintf(file," %d",  bin[i]);                        fprintf(file,"\n");
 fprintf(file,"#! MIN");  for(i=0;i<grid->ncv;i++) fprintf(file," %lf", grid->min[i]);                  fprintf(file,"\n"); 
 fprintf(file,"#! MAX");  for(i=0;i<grid->ncv;i++) fprintf(file," %lf", max[i]);                        fprintf(file,"\n");
 fprintf(file,"#! PBC");  for(i=0;i<grid->ncv;i++) fprintf(file," %d",  grid->period[i]);               fprintf(file,"\n");

// GRID
 for(i=0;i<grid->size;i++){
  for(j=0;j<grid->ncv;j++) {
   xx[j] = grid->min[j] + grid->dx[j] * grid->one2multi_full[i][j];  
   fprintf(file," %lf ",xx[j]);
  }
  fprintf(file," %lf ",grid->pot[i]/mtd_data.eunit);
  for(j=0;j<grid->ncv;j++) fprintf(file," %lf ",grid->force[i][j]/mtd_data.eunit); fprintf(file,"\n");
  if(grid->one2multi_full[i][0]==(grid->bin[0]-1)) fprintf(file,"\n");
 }

// Deallocation
 free_1dr_array_alloc(xx);
 free_1dr_array_alloc(max);
 free_1di_array_alloc(index_nd);  
 free_1di_array_alloc(bin);

// Final stuff
 fflush(file);
 fclose(file);

}
//-------------------------------------------------------------------------------------------
// read a GRID from file
void  PREFIX grid_read_fromfile(struct grid_s *grid, int bias)
{
 struct grid_s tmpgrid;
 int      i, j, with_force=0, isave, line, m;
 int      header, *index_nd, index_1d;
 char     str1[30], str2[30], str3[30];
 char     *str, stringa[800];
 FILE     *file=NULL;
 real   *ff, *xx, Vgrid=0., Vp, Vm;
 double  tmp;

// open grid file for reading
 file = fopen(grid->r_file, "r");
 if(!file) plumed_error("Cannot read GRID file!\n");

// read HEADER
 header = 0;
 fprintf(mtd_data.fplog,"** READING GRID FROM FILE %s\n",grid->r_file);
// first you need NVAR for allocation
 str = fgets(stringa, 800, file);
 while(1){
   if(sscanf(str,"%s %s%n",str1,str2,&m)>0){
    if(strcmp(str1,"#!") == 0){
     str +=m;
     if(strcmp(str2,"NVAR") == 0)  {sscanf(str,"%d",&(tmpgrid.ncv)); header +=1;}
    } else if(strcmp(str1,"#") != 0) break;
   }
   str = fgets(stringa, 800, file);
 }
 // checking for missing or wrong data 
 if(header!=1) plumed_error("Missing or wrong data in GRID header!\n");
// and then parse the rest
 rewind(file);
 str = fgets(stringa, 800, file);
 while(1){
   if(sscanf(str,"%s %s%n",str1,str2,&m)>0){
    if(strcmp(str1,"#!") == 0){
     str +=m;
     if(strcmp(str2,"FORCE") == 0) {sscanf(str,"%s%n",str3,&m); with_force=atoi(str3); header +=1;str +=m;}
     if(strcmp(str2,"TYPE") == 0)  {for(i=0;i<tmpgrid.ncv;i++) {sscanf(str,"%s%n",str3,&m); tmpgrid.index[i]  = atoi(str3);str +=m;} header +=1;}
     if(strcmp(str2,"BIN") == 0)   {for(i=0;i<tmpgrid.ncv;i++) {sscanf(str,"%s%n",str3,&m); tmpgrid.bin[i]    = atoi(str3);str +=m;} header +=1;}
     if(strcmp(str2,"MIN") == 0)   {for(i=0;i<tmpgrid.ncv;i++) {sscanf(str,"%s%n",str3,&m); tmpgrid.min[i]    = atof(str3);str +=m;} header +=1;}
     if(strcmp(str2,"MAX") == 0)   {for(i=0;i<tmpgrid.ncv;i++) {sscanf(str,"%s%n",str3,&m); tmpgrid.max[i]    = atof(str3);str +=m;} header +=1;}
     if(strcmp(str2,"PBC") == 0)   {for(i=0;i<tmpgrid.ncv;i++) {sscanf(str,"%s%n",str3,&m); tmpgrid.period[i] = atoi(str3);str +=m;} header +=1;}
    } else if(strcmp(str1,"#") == 0) printf("   COMMENT: %s",stringa);
    else break;
   }
   str = fgets(stringa, 800, file);
 }

// checking for missing or wrong data 
 if(header!=7) plumed_error("Missing or wrong data in GRID header!\n");
// compare with  grid.ncv
 if(grid->ncv!=tmpgrid.ncv) plumed_error("Inconsistency between NVAR on file and in the PLUMED input file\n");
// compare with grid.index
 for(i=0;i<tmpgrid.ncv;i++) if(tmpgrid.index[i]!=(colvar.type_s[grid->index[i]])) 
  plumed_error("Inconsistency between CV TYPES on file and in the PLUMED input file\n"); 
// printout HEADER 
 fprintf(mtd_data.fplog,"   NVAR :: %d \n",tmpgrid.ncv);
 fprintf(mtd_data.fplog,"   BIN  :: "); for(i=0;i<tmpgrid.ncv;i++) fprintf(mtd_data.fplog,"%d ",tmpgrid.bin[i]); fprintf(mtd_data.fplog,"\n");
 fprintf(mtd_data.fplog,"   MIN  :: "); for(i=0;i<tmpgrid.ncv;i++) fprintf(mtd_data.fplog,"%lf ",tmpgrid.min[i]); fprintf(mtd_data.fplog,"\n");
 fprintf(mtd_data.fplog,"   MAX  :: "); for(i=0;i<tmpgrid.ncv;i++) fprintf(mtd_data.fplog,"%lf ",tmpgrid.max[i]); fprintf(mtd_data.fplog,"\n");
 fprintf(mtd_data.fplog,"   PBC  :: "); for(i=0;i<tmpgrid.ncv;i++) if(tmpgrid.period[i]==0) fprintf(mtd_data.fplog,"OFF "); else fprintf(mtd_data.fplog,"ON "); 
 fprintf(mtd_data.fplog,"\n");

// if reading an external potential we can start initialize something
// No need if bias==1 since grid is initialized in read_restraint
 if(bias==0){                           
  for(i=0;i<grid->ncv;i++){
   grid->bin[i]    = tmpgrid.bin[i];
   grid->min[i]    = tmpgrid.min[i];
   grid->max[i]    = tmpgrid.max[i];
   grid->period[i] = tmpgrid.period[i];
  }
  grid_initialize(grid);
 }

// copying grid.index into tmpgrid.index
 for(i=0;i<tmpgrid.ncv;i++) tmpgrid.index[i]=grid->index[i];
// and initializing tmpgrid
 grid_initialize(&tmpgrid);

// allocating temp arrays
 xx       = float_1d_array_alloc(tmpgrid.ncv);
 ff       = float_1d_array_alloc(tmpgrid.ncv);
 index_nd = int_1d_array_alloc(tmpgrid.ncv);

// now parsing the grid for potential and forces 
 line = 0;
 while(1){                                     
   if(sscanf(str,"%s %s",str1,str2)>0){
     j=0;
     str =(char *) strtok(stringa," \t");
     while (str != NULL)
     { 
      if(j>=0 && j<tmpgrid.ncv) { sscanf(str, "%lf", &tmp); xx[j]=(real)tmp;}
      if(j==tmpgrid.ncv)       { sscanf(str, "%lf", &tmp); Vgrid = (real)tmp;}
      if(j>tmpgrid.ncv && j<=2*tmpgrid.ncv && with_force==1) {sscanf(str, "%lf", &tmp); ff[j-tmpgrid.ncv-1]=(real)tmp;}
      str =(char *) strtok(NULL, " \t");
      j++;
     }
// find multi dimensional index
     for(i=0;i<tmpgrid.ncv;i++) index_nd[i] = floor((xx[i]+tmpgrid.dx[i]/2.-tmpgrid.min[i])/tmpgrid.dx[i]);
// and mono-dimensional
     index_1d=grid_multi2one(&tmpgrid,index_nd);
     if(index_1d!=line) plumed_warn("GRID on file is not in the usual PLUMED format");
// Storing potential...
     tmpgrid.pot[index_1d]=Vgrid*mtd_data.eunit;
// ...and forces 
     if(with_force==1) for(i=0;i<tmpgrid.ncv;i++) tmpgrid.force[index_1d][i]=ff[i]*mtd_data.eunit; 
// new line
     line++;
   }
     str = fgets(stringa, 800, file);
     if(str == NULL) break;
 }

// check total size and line
  if(line!=tmpgrid.size) plumed_error("GRID entries on file are not consistent with the declared dimension \n");

// if derivatives are missing, finite differences... 
  if(with_force==1) fprintf(mtd_data.fplog,"** FORCE DATA ARE PRESENT ON FILE\n");
  else {
   fprintf(mtd_data.fplog,"** NO FORCE DATA ON FILE: FINITE DIFFERENCES\n");
   for(i=0;i<tmpgrid.size;i++){
    for(j=0;j<tmpgrid.ncv;j++) index_nd[j] = tmpgrid.one2multi_full[i][j];
    for(j=0;j<tmpgrid.ncv;j++){
       isave=index_nd[j]; index_nd[j] += 1;
       if(index_nd[j]==tmpgrid.bin[j]  && tmpgrid.period[j]==0) {tmpgrid.force[i][j]=0.0; continue;}
       if(index_nd[j]==tmpgrid.bin[j]  && tmpgrid.period[j]==1) index_nd[j] = 0; 
       index_1d=grid_multi2one(&tmpgrid,index_nd); 
       Vp=tmpgrid.pot[index_1d];
       index_nd[j]=isave; index_nd[j] -= 1;
       if(index_nd[j]==-1  && tmpgrid.period[j]==0) {tmpgrid.force[i][j]=0.0; continue;}
       if(index_nd[j]==-1  && tmpgrid.period[j]==1) index_nd[j] = tmpgrid.bin[j]-1;
       index_1d=grid_multi2one(&tmpgrid,index_nd);
       Vm=tmpgrid.pot[index_1d]; 
       index_nd[j]=isave;
       tmpgrid.force[i][j]=-(Vp-Vm)/2.0/tmpgrid.dx[j];
    }
   }
  }

// Now tmpgrid is complete. Time to clone it to grid.
// And interpolate if necessary. 
 grid_clone(&tmpgrid, grid);

// Deallocation
 free_1dr_array_alloc(ff);
 free_1dr_array_alloc(xx);
 free_1di_array_alloc(index_nd);
 free_2di_array_alloc(tmpgrid.one2multi_full,tmpgrid.size);
 free_1dr_array_alloc(tmpgrid.pot);
 free_2dr_array_alloc(tmpgrid.force, tmpgrid.size);

// Final stuff
 fclose(file);
 fprintf(mtd_data.fplog,"\n");
} 
//-------------------------------------------------------------------------------------------
void PREFIX grid_clone(struct grid_s *grid1, struct grid_s *grid2)
{

 int      i, j, just_copy, out_grid;
 real   xx[nconst_max], ff[nconst_max]; 

// A copy is enough ??
// - check number of bins
// - check boundaries
 just_copy=1;
 for(i=0;i<grid1->ncv;i++) {
  if(grid1->bin[i]!=grid2->bin[i]) just_copy=0;
  if(fabs(grid1->min[i]-grid2->min[i])>0.00001) just_copy=0;
  if(fabs(grid1->max[i]-grid2->max[i])>0.00001) just_copy=0;
 }

 fprintf(mtd_data.fplog,"** CLONING GRID "); if(just_copy==0) fprintf(mtd_data.fplog,"AND INTERPOLATING"); fprintf(mtd_data.fplog,"\n");

 if(just_copy==1) for(i=0;i<grid1->size;i++) {
   grid2->pot[i]=grid1->pot[i]; 
   for(j=0;j<grid1->ncv;j++) grid2->force[i][j]=grid1->force[i][j];} 
 else { // Need interpolation
  for(i=0;i<grid2->size;i++){
   out_grid = 0;
   for(j=0;j<grid2->ncv;j++) {
    ff[grid2->index[j]]=0.0;
    xx[grid2->index[j]] = grid2->min[j] + grid2->dx[j] * grid2->one2multi_full[i][j];
    if(xx[grid2->index[j]]<grid1->min[j] || xx[grid2->index[j]]>=grid1->max[j]-grid1->dx[j]) out_grid = 1;
   } 
   if(out_grid==1) {grid2->pot[i]=0.0; for(j=0;j<grid2->ncv;j++) grid2->force[i][j]=0.0;} 
   else{
    grid2->pot[i]=grid_getstuff(grid1,xx,ff);
    for(j=0;j<grid2->ncv;j++) grid2->force[i][j]=ff[grid2->index[j]];
   }
  }
 }

}
//-------------------------------------------------------------------------------------------
// Interpolation with a (sort of) cubic spline.
// The function is built as a sum over the nearest neighbours (i.e. 2 in 1d, 4 in 2d, 8 in 3d,...).
// Each neighbour contributes with a polynomial function which is a product of single-dimensional polynomials,
// written as functions of the distance to the neighbour in units of grid spacing
// Each polynomial is proportional to:
// (1-3x^2+2x^3)  + Q (x-2x^2+x^3)
// * its value and derivative in +1 are zero
// * its value in 0 is 1
// * its derivative in 0 is Q
// so, Q is chosen as the desired derivative at the grid point divided by the value at the grid point
// and the final function is multiplied times the value at the grid point.
//
// It works perfectly, except when the tabulated function is zero (there is a special case).
// Maybe one day I will learn the proper way to do splines...
// Giovanni

real PREFIX spline(int ndim,real *dx,real *where,real *tabf,real *tabder,int* stride,real *der){
// ndim:   dimensionality
// dx:     delta between grid points
// where:  location relative to the floor grid point (always between 0 and dx)
// tabf:   table with function, already pointed at the floor grid point
// tabder: table with minus gradients (the fastest running index is the dimension index), already pointed at the floor grid point
// stride: strides to the next point on the tabf array.
//         note that, in case of PBC, this stride should corrispond to a backward jump of (N-1) points,
//         where N is the number of points in the domain. 
//         also note that the corrisponding strides for tabder can be obtained multipling times ndim
// der:    in output, the minus gradient.

  int idim;
  int npoints,ipoint;
  real X;
  real X2;
  real X3;
  int x0[nconst_max];;
  real fd[nconst_max];
  real C[nconst_max];
  real D[nconst_max];
  int  tmp,shift;
  real f;

  npoints=1; for(idim=0;idim<ndim;idim++) npoints*=2; // npoints=2**ndim

// reset
  f=0;
  for(idim=0;idim<ndim;idim++) der[idim]=0;

// loop over neighbour points:
  for(ipoint=0;ipoint<npoints;ipoint++){

// find coordinate of neighbour point (x0) and shift
    tmp=ipoint;
    shift=0;
    for(idim=0;idim<ndim;idim++){
      x0[idim]=tmp%2; tmp/=2;
      shift+=stride[idim]*x0[idim];
    }
//fprintf(stderr,"%i\n",shift);

// reset contribution from this point:
    real ff;
    ff=1.0;

    for(idim=0;idim<ndim;idim++){
      X=fabs(where[idim]/dx[idim]-x0[idim]);
      X2=X*X;
      X3=X2*X;
      real yy;
      if(fabs(tabf[shift])<0.0000001) yy=0.0;
      else yy=tabder[shift*ndim+idim]/tabf[shift];
                                       // il - e per -derivata
      C[idim]=(1-3*X2+2*X3) - (x0[idim]?-1:1)*yy*(X-2*X2+X3)*dx[idim];
      D[idim]=( -6*X +6*X2) - (x0[idim]?-1:1)*yy*(1-4*X +3*X2)*dx[idim]; // d / dX
      D[idim]*=(x0[idim]?-1:1)/dx[idim]; // chain rule (to where)
      ff*=C[idim];
    }
    for(idim=0;idim<ndim;idim++) {
      int idim1;
      fd[idim]=D[idim];
      for(idim1=0;idim1<ndim;idim1++) if(idim1!=idim) fd[idim]*=C[idim1];
    }
    
    f+=tabf[shift]*ff;
    for(idim=0;idim<ndim;idim++) der[idim]+=tabf[shift]*fd[idim];
  }
  return f;
};

// JFD>
//-------------------------------------------------------------------------------------------
// Transition-tempered metadynamics implementation details.
// This section provides all the utilities necessary to find the
// current best guess of the simulation's bias at a transition state
// between two wells and the function that computes that best guess. 

// There are four
// primary utilities provided. First is a lattice data structure that is used to
// interpret the grid data structure that PLUMED stores the bias in. Second is a 
// data structure for functions defined on that lattice, like the bias. Third is a
// max-heap data structure that is used as a priority queue in finding the maximum of 
// the minimal bias along any paths between two points. The final utility will find
// maxima of minimal paths for a provided positive function defined on the grid points. 

// The first fully-defined function calculates the transition bias for the 
// hills engine. This is the only function that uses PLUMED functions other
// than plumed_error.

double PREFIX transition_bias_ND() {
    double transition_bias, least_transition_bias;
    lattice *lat;
    lattice_array *lat_pot;
    size_t source, sink;
    size_t *d_sizes;
    size_t *d_bcs;
    size_t dims;
    size_t i, j;

    // Use the grid and bias to define convenient minimal data structures.
    dims = bias_grid.ncv;
    d_sizes = (size_t *)malloc(dims * sizeof(size_t));
    d_bcs = (size_t *)malloc(dims * sizeof(size_t));
    for (i = 0; i < dims; i++) {
        d_sizes[i] = bias_grid.bin[i];
        if (bias_grid.period[i]) {
            d_bcs[i] = TTMETAD_LATTICE_PBC;
        } else if (!bias_grid.period[i]) {
            d_bcs[i] = TTMETAD_LATTICE_RBC;
        }
    }

    lat = make_lattice(dims, d_sizes, d_bcs);
    lat_pot = make_const_lattice_array(lat);
    set_const_lattice_array(lat_pot, bias_grid.pot);


    if (!colvar.transition_wells_converted) {
        int *index_nd;

        // Allocate an array for the multidimensional grid index of each well.
        index_nd = int_1d_array_alloc(dims);
        for (i = 0; i < colvar.n_transition_wells; i++) {

            // Convert to multidimensional grid indices.
            for(j = 0; j < dims; j++) {
                if((colvar.transition_wells[i][j]<bias_grid.min[j] || colvar.transition_wells[i][j]>=bias_grid.max[j]) && !bias_grid.period[j]) plumed_error("Transition well is outside the GRID!. Please increase GRID size.");
                if(bias_grid.period[j])  colvar.transition_wells[i][j] -= bias_grid.lbox[j] * rint(colvar.transition_wells[i][j]/bias_grid.lbox[j]);
                index_nd[j] = floor((colvar.transition_wells[i][j]-bias_grid.min[j])/bias_grid.dx[j]);
            }

            // Convert to one-dimensional index and record it.
            colvar.transition_well_indices[i] = grid_multi2one(&bias_grid, index_nd);
        }
        colvar.transition_wells_converted = 1;
        free(index_nd);
    }

    // Find the transition bias using those variables.
    // The transition is the maximum of minimal maxima of paths between pairs of wells.
    least_transition_bias = find_maximal_path_minimum(lat_pot, colvar.transition_well_indices[1], colvar.transition_well_indices[0]);
    for (i = 0; i < colvar.n_transition_wells; i++) {
        for (j = 0; j < i; j++) {
            if (i == 1 && j == 0) continue;    // Already calculated.
            source = colvar.transition_well_indices[i];
            sink = colvar.transition_well_indices[j];
            transition_bias = find_maximal_path_minimum(lat_pot, source, sink);
            least_transition_bias = fmin(transition_bias, least_transition_bias);
        }
    }
    
    // Clean up.
    free_const_lattice_array(lat_pot);
    free_lattice(lat);
    free(d_bcs);
    free(d_sizes);
    
    // Threshold the return value and send it on.
    transition_bias = fmax(0.0, least_transition_bias - colvar.ttthreshold);
    return transition_bias;
}

PREFIX lattice *PREFIX make_lattice(size_t dims, size_t *d_sizes, size_t *d_bcs) {
    lattice *new_lat;
    size_t i;
    size_t stride = 1;

    new_lat = (lattice *)malloc(sizeof(lattice));

    new_lat->dimensions = dims;

    new_lat->dim_sizes = (size_t *)malloc(dims * sizeof(size_t));
    new_lat->dim_strides = (size_t *)malloc(dims * sizeof(size_t));
    new_lat->dim_bcs = (size_t *)malloc(dims * sizeof(size_t));

    /* Set the size of the lattice in each dimension */
    for (i = 0; i < dims; ++i) {
        new_lat->dim_sizes[i] = d_sizes[i];
    }

    /* Set the lattice strides and total size for the lattice */
    for (i = 0; i < dims; ++i) {
        new_lat->dim_strides[i] = stride;
        stride *= new_lat->dim_sizes[i];
    }
    new_lat->total_size = stride;

    /* Set the boundary conditions of the lattice in each dimension */
    for (i = 0; i < dims; ++i) {
        new_lat->dim_bcs[i] = d_bcs[i];
    }

    return new_lat;
}

void PREFIX free_lattice(PREFIX lattice *lat) {
    free(lat->dim_sizes);
    free(lat->dim_strides);
    free(lat->dim_bcs);
    free(lat);
}

size_t *PREFIX one_d_index_to_multi_d(PREFIX lattice *lat, size_t index) {
    
    size_t i, last_dimension;
    size_t *multi_d_index;

    multi_d_index = (size_t *)malloc(lat->dimensions * sizeof(size_t));
    last_dimension = lat->dimensions - 1;
    for (i = 0; i < last_dimension; ++i) {
        multi_d_index[i] = (index / lat->dim_strides[i]) % lat->dim_sizes[i + 1];
    }
    multi_d_index[last_dimension] = (index / lat->dim_strides[last_dimension]);
    return multi_d_index;
}

size_t PREFIX multi_d_index_to_one_d(PREFIX lattice *lat, size_t *index) {
    
    size_t i;
    size_t single_index = 0;
    
    for (i = 0; i < lat->dimensions; ++i) {
        single_index += index[i] * lat->dim_strides[i];
    }
    
    return single_index;
}

void PREFIX find_neighbors(PREFIX lattice *lat, size_t index, size_t *n_neighs, size_t *neighs) {
    
    size_t i;
    size_t *dim_indices;
    dim_indices = one_d_index_to_multi_d(lat, index);
    if (neighs == NULL) {
        exit(EXIT_FAILURE);
    }
    *n_neighs = 0;

    for (i = 0; i < lat->dimensions; ++i) {
        if (dim_indices[i] == 0) {
            if (lat->dim_bcs[i] == TTMETAD_LATTICE_PBC) {
                dim_indices[i] = 1;
                neighs[*n_neighs] = multi_d_index_to_one_d(lat, dim_indices);
                *n_neighs += 1;
                dim_indices[i] = lat->dim_sizes[i] - 1;
                neighs[*n_neighs] = multi_d_index_to_one_d(lat, dim_indices);
                *n_neighs += 1;
                dim_indices[i] = 0;
            } else if (lat->dim_bcs[i] == TTMETAD_LATTICE_RBC) {
                dim_indices[i] = 1;
                neighs[*n_neighs] = multi_d_index_to_one_d(lat, dim_indices);
                *n_neighs += 1;
                dim_indices[i] = 0;
            }
        } else if (dim_indices[i] == (lat->dim_sizes[i] - 1)) {
            if (lat->dim_bcs[i] == TTMETAD_LATTICE_PBC) {
                dim_indices[i] = 0;
                neighs[*n_neighs] = multi_d_index_to_one_d(lat, dim_indices);
                *n_neighs += 1;
                dim_indices[i] = lat->dim_sizes[i] - 2;
                neighs[*n_neighs] = multi_d_index_to_one_d(lat, dim_indices);
                *n_neighs += 1;
                dim_indices[i] = lat->dim_sizes[i] - 1;
            } else if (lat->dim_bcs[i] == TTMETAD_LATTICE_RBC) {
                dim_indices[i] = lat->dim_sizes[i] - 2;
                neighs[*n_neighs] = multi_d_index_to_one_d(lat, dim_indices);
                *n_neighs += 1;
                dim_indices[i] = lat->dim_sizes[i] - 1;
            }
        } else {
            dim_indices[i] += 1;
            neighs[*n_neighs] = multi_d_index_to_one_d(lat, dim_indices);
            *n_neighs += 1;
            dim_indices[i] -= 2;
            neighs[*n_neighs] = multi_d_index_to_one_d(lat, dim_indices);
            *n_neighs += 1;
            dim_indices[i] += 1;
        }
    }
    free(dim_indices);
}


/* Lattice-dependent data structure routines */

PREFIX lattice_array *PREFIX make_lattice_array(PREFIX lattice *lat) {
    
    lattice_array *new_lat_arr;
    size_t i;
    
    /* Allocate the struct. */
    new_lat_arr = (lattice_array *)malloc(sizeof(lattice_array));
    
    /* Set which lattice pointer to use. */
    new_lat_arr->index_lat = lat;
    
    /* Allocate and initialize the array to zero. */
    new_lat_arr->vals = (double *)malloc(lat->total_size * sizeof(double));
    for (i = 0; i < lat->total_size; ++i) {
        new_lat_arr->vals[i] = 0;
    }
    
    return new_lat_arr;
}

PREFIX lattice_array *PREFIX make_const_lattice_array(PREFIX lattice *lat) {
    
    lattice_array *new_lat_arr;
    
    /* Allocate the struct. */
    new_lat_arr = (lattice_array *)malloc(sizeof(lattice_array));
    
    /* Set which lattice pointer to use. */
    new_lat_arr->index_lat = lat;
    
    /* Do not allocate or initialize the array to zero. */
    new_lat_arr->vals = NULL;
    
    return new_lat_arr;
}

void PREFIX free_lattice_array(PREFIX lattice_array *lat_arr) {
    free(lat_arr->vals);
    free(lat_arr);
}

void PREFIX free_const_lattice_array(PREFIX lattice_array *lat_arr) {
    free(lat_arr);
}


void PREFIX set_lattice_array(PREFIX lattice_array *lat_arr, double *vals) {
    size_t i;
    for (i = 0; i < lat_arr->index_lat->total_size; ++i) {
        lat_arr->vals[i] = vals[i];
    }
}

void PREFIX set_const_lattice_array(PREFIX lattice_array *lat_arr, double *vals) {
    lat_arr->vals = vals;
}


/* Min heap routines for indexed values */

PREFIX indexed_value_max_heap *PREFIX make_indexed_value_max_heap(size_t max_size) {

    indexed_value_max_heap *new_ival_mheap;
    size_t i;

    new_ival_mheap = (indexed_value_max_heap *)malloc(sizeof(indexed_value_max_heap));
    new_ival_mheap->max_size = max_size;
    new_ival_mheap->curr_size = 0;
    new_ival_mheap->indices = (size_t *)malloc(max_size * sizeof(size_t));
    for (i = 0; i < max_size; ++i) {
        new_ival_mheap->indices[i] = 0;
    }
    new_ival_mheap->vals = (double *)malloc(max_size * sizeof(double));
    for (i = 0; i < max_size; ++i) {
        new_ival_mheap->vals[i] = 0.0;
    }
    return new_ival_mheap;
}

void PREFIX free_indexed_value_max_heap(PREFIX indexed_value_max_heap *ival_mheap) {
    free(ival_mheap->indices);
    free(ival_mheap->vals);
    free(ival_mheap);
}

double PREFIX get_max(PREFIX indexed_value_max_heap *ival_mheap) {
    return ival_mheap->vals[0];
}

size_t PREFIX get_max_index(PREFIX indexed_value_max_heap *ival_mheap) {
    return ival_mheap->indices[0];
}

void PREFIX swap_heap_entries(PREFIX indexed_value_max_heap *ival_mheap, size_t index_one, size_t index_two) {
    
    size_t temp_index;
    double temp_val;
    
    temp_val = ival_mheap->vals[index_one];
    temp_index =  ival_mheap->indices[index_one];
    ival_mheap->vals[index_one] = ival_mheap->vals[index_two];
    ival_mheap->indices[index_one] = ival_mheap->indices[index_two];
    ival_mheap->vals[index_two] = temp_val;
    ival_mheap->indices[index_two] = temp_index;
}

void PREFIX insert_indexed_value(PREFIX indexed_value_max_heap *ival_mheap, size_t index, double value) {
    
    size_t last_index, pivot_index, parent_index;
    int improperly_heaped;

    /* Insert the data at the end of the heap. */
    if (ival_mheap->curr_size >= ival_mheap->max_size) {
        char buf[1024];
        sprintf(buf,"Error finding transition bias: heap of size %ld is already full of %ld entries.\n", ival_mheap->max_size, ival_mheap->curr_size);
        plumed_error(buf);
    }
    last_index = ival_mheap->curr_size;
    ival_mheap->indices[last_index] = index;
    ival_mheap->vals[last_index] = value;
    ival_mheap->curr_size++;
    improperly_heaped = 1;

    /* Re-structure the new state to get a proper heap. */
    pivot_index = last_index;
    while (improperly_heaped) {
        parent_index = pivot_index / 2;
        if (ival_mheap->vals[pivot_index] > ival_mheap->vals[parent_index]) {
            swap_heap_entries(ival_mheap, pivot_index, parent_index);
            pivot_index = parent_index;
        } else {
            improperly_heaped = 0;
        }
    }
}

void PREFIX delete_max(PREFIX indexed_value_max_heap *ival_mheap) {
    
    size_t last_index, pivot_index, child_index_1, child_index_2;
    int improperly_heaped;
    
    /* Overwrite the current minimum with the last entry in the heap
       and zero the last entry. */
    if (ival_mheap->curr_size > 1) {
        last_index = ival_mheap->curr_size - 1;
        ival_mheap->indices[0] = ival_mheap->indices[last_index];
        ival_mheap->vals[0] = ival_mheap->vals[last_index];
        ival_mheap->indices[last_index] = 0;
        ival_mheap->vals[last_index] = 0.0;
        ival_mheap->curr_size--;
        improperly_heaped = 1;

        /* Re-structure the new state to get a proper heap. */
        pivot_index = 0;
        while (improperly_heaped) {
            /* Examine the child elements of the pivot element. */
            child_index_1 = 2 * pivot_index + 1;
            child_index_2 = 2 * pivot_index + 2;
        
            /* If the element has no children, the heap is ready. */
            if (child_index_1 >= ival_mheap->curr_size) {
                improperly_heaped = 0;
        
            /* If the element has only one child, consider swapping them. */
            } else if (child_index_2 >= ival_mheap->curr_size) {
                if (ival_mheap->vals[child_index_1] > ival_mheap->vals[pivot_index]) {
                    swap_heap_entries(ival_mheap, pivot_index, child_index_1);
                    pivot_index = child_index_1;
                }
                improperly_heaped = 0;

            /* If it has two children, swap the pivot element with its greatest 
               child only if they are mis-ordered. */
            } else {
                if (ival_mheap->vals[pivot_index] >= ival_mheap->vals[child_index_1] &&
                    ival_mheap->vals[pivot_index] >= ival_mheap->vals[child_index_2]) {
                    improperly_heaped = 0;
                } else {
                    if (ival_mheap->vals[child_index_1] > ival_mheap->vals[child_index_2]) {
                        swap_heap_entries(ival_mheap, pivot_index, child_index_1);
                        pivot_index = child_index_1;
                    } else {
                        swap_heap_entries(ival_mheap, pivot_index, child_index_2);
                        pivot_index = child_index_2;
                    }
                }
            }
        }
    } else if (ival_mheap->curr_size == 1) {
        ival_mheap->indices[0] = 0;
        ival_mheap->vals[0] = 0.0;
        ival_mheap->curr_size--;
    } else {
        plumed_error("Error finding transition bias: deleting from heap that is already empty.");
    }
}

double PREFIX find_maximal_path_minimum(lattice_array *lat_arr, size_t source, size_t sink) {
    
    int path_not_found = 1;
    double maximal_minimum = 0.0;
    size_t i;
    lattice *lat = lat_arr->index_lat;
    size_t lattice_size = lat->total_size;
    double *mins_to_each;
    size_t *steps_back;
    indexed_value_max_heap *next_steps;

    size_t curr_point, n_neighs, n_examined;
    size_t *neighs;
    size_t max_path_length = 0;

    FILE *debugfp;
    clock_t path_search_start_time, path_elapsed_search_time;
    int msecs;

    if (source >= lat->total_size) {
        char buf[1024];
        sprintf(buf,"Error finding transition bias: source %ld is not in the bias domain of size %ld.\n", source, lat->total_size);
        plumed_error(buf);
    }
    if (source >= lat->total_size) {
        char buf[1024];
        sprintf(buf,"Error finding transition bias: sink %ld is not in the bias domain of size %ld.", sink, lat->total_size);
        plumed_error(buf);
    }

    mins_to_each = (double *)malloc(lattice_size * sizeof(double));
    for (i = 0; i < lattice_size; i++) {
        mins_to_each[i] = -1.0;
    }
    steps_back = (size_t *)malloc(lattice_size * sizeof(size_t));
    for (i = 0; i < lattice_size; i++) {
        steps_back[i] = lattice_size;
    }
    next_steps = make_indexed_value_max_heap(lattice_size);
    neighs = (size_t *)malloc( 2 * lat->dimensions * sizeof(size_t));
    for (i = 0; i < lat->dimensions; ++i) {
        max_path_length += lat->dim_sizes[i];
    }
    
    /* 
     * If debugging, open the debug file to append new data and 
     * write the start of a line to the debug file. Also record
     * the time the path search began.
     *
     */
    if (logical.ttdebug) {
      debugfp = fopen(colvar.ttdebug_file, "a");
      fprintf(debugfp, "Searching for path from %ld to %ld.\n", source, sink);
      path_search_start_time = clock();
    }

    /* Use Dijkstra's algorithm using the minimum path value as the 
       total inverse cost of the path. */
    mins_to_each[source] = lat_arr->vals[source];
    insert_indexed_value(next_steps, source, lat_arr->vals[source]);
    n_examined = 1;
    
    while (path_not_found) {
        
        /* Get the current best-looking path step to take. */
        curr_point = get_max_index(next_steps);
        delete_max(next_steps);
        
        /* If that takes us to the sink, we're done. */
        if (curr_point == sink) {
            path_not_found = 0;
            maximal_minimum = mins_to_each[sink];
            break;
        /* Also, if the current maximal minimum is the lowest possible minimum we're done. */
        } else if (mins_to_each[curr_point] == 0.0) {
            maximal_minimum = 0.0;
            break;
        }
        
        /* Find its neighbors. */
        find_neighbors(lat, curr_point, &n_neighs, neighs);
        /* Add the neighbors (with relevant path minima) to the step queue */ 
        for (i = 0; i < n_neighs; ++i) {
            if (mins_to_each[neighs[i]] == -1.0) {
                mins_to_each[neighs[i]] = fmin(mins_to_each[curr_point], lat_arr->vals[neighs[i]]);
                steps_back[neighs[i]] = curr_point;
                insert_indexed_value(next_steps, neighs[i], mins_to_each[neighs[i]]);
            }
        }
        n_examined++;
        if (n_examined >= lattice_size) {
            plumed_error("Faulty loop; the path-finding algorithm is retracing its steps.\n");
        }
        /* Move on to the next best-looking next step along any path. */
    }

    /* If debugging, write the conclusion to the debug file. */
    if (logical.ttdebug) {
      path_elapsed_search_time = (clock() - path_search_start_time);
      msecs = path_elapsed_search_time * 1000 / CLOCKS_PER_SEC;
      fprintf(debugfp, "Path search complete in %ld steps, taking %d msec.\n", n_examined, msecs);
      fclose(debugfp);
    }

    free(mins_to_each);
    free(steps_back);
    free_indexed_value_max_heap(next_steps);
    free(neighs);

    return maximal_minimum;
}

// <JFD
