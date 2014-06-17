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

/*
 * This is an implementation of experiment directed simulations.
 *
 */



/*
 * set-up the parameters. hard cou
 */ 
void PREFIX eds_init(int cv_number, real update_period, 
		       real simtemp, int seed,
		       int b_hard_coupling_range, 
		       int* cv_map,
		       t_eds* eds) {

  eds->centers = (real*) calloc(cv_number, sizeof(real));
  eds->means = (real*) calloc(cv_number, sizeof(real));
  eds->ssd = (real*) calloc(cv_number, sizeof(real));
  eds->max_coupling_range = (real*) calloc(cv_number, sizeof(real));
  eds->max_coupling_rate = (real*) calloc(cv_number, sizeof(real));
  eds->set_coupling = (real*) calloc(cv_number, sizeof(real));
  eds->current_coupling = (real*) calloc(cv_number, sizeof(real));
  eds->coupling_rate = (real*) calloc(cv_number, sizeof(real));
  eds->coupling_accum = (real*) calloc(cv_number, sizeof(real));

  eds->cv_number = cv_number;
  
  eds->cv_map = cv_map;  
  eds->simtemp = simtemp;  
  eds->seed = seed;

  eds->update_period = update_period;
  eds->update_calls = 0;
  eds->b_equilibration = 1;
  eds->b_hard_coupling_range = b_hard_coupling_range;
 
}

void PREFIX eds_free(t_eds* eds) {

  free(eds->centers);
  free(eds->means);
  free(eds->ssd);
  free(eds->max_coupling_rate);
  free(eds->max_coupling_range);
  free(eds->set_coupling);
  free(eds->current_coupling);
  free(eds->coupling_rate);
  free(eds->coupling_accum);

}

/*
 * Run update of calculation and apply forces
 */
real PREFIX eds_engine(real* ss0, real* force, 
		       t_eds* eds, real boltz) {

  real bias_energy = 0.0;
  int i;  

  eds->update_calls++;
  int b_finished_equil_flag = 1;
  real delta;
    
  //apply forces for this setp and calculat energies
  for(i = 0; i < eds->cv_number; i++) {
    force[eds->cv_map[i]] = eds->current_coupling[i] / eds->centers[i] * ss0[eds->cv_map[i]];
    bias_energy += eds->current_coupling[i] / eds->centers[i] * (ss0[eds->cv_map[i]] - eds->centers[i]);

    //are we updating the bias?
    if(eds->update_period == 0)
      continue;
    
    //if we aren't waiting for the bias to equilibrate
    if(!eds->b_equilibration) {

      //Welford, West, and Hanso online variance method
      delta = ss0[eds->cv_map[i]] - eds->means[i];
      eds->means[i] += delta / eds->update_calls;
      eds->ssd[i] += delta * (ss0[eds->cv_map[i]] - eds->means[i]);
    } else {
      //equilibrating
      //check if we've reached the setpoint
      if(eds->coupling_rate[i] == 0 || pow(eds->current_coupling[i] - eds->set_coupling[i], 2) < pow(eds->coupling_rate[i], 2)) {
	b_finished_equil_flag &= 1;
      }
      else {
	eds->current_coupling[i] += eds->coupling_rate[i];
	b_finished_equil_flag = 0;
      }
    }
    
    //still in for loop. Now we update our max coupling range
    //if we're allowing that
    if(!eds->b_hard_coupling_range && fabs(eds->current_coupling[i]) > 
       eds->max_coupling_range[i]) {
      eds->max_coupling_range[i] *= 1.25;
    }        
  }

  //reduce all the flags 
  if(eds->b_equilibration && b_finished_equil_flag) {
    eds->b_equilibration = false;
    eds->update_calls = 0;
  }


  //Now we update coupling constant, if necessary
  if(!eds->b_equilibration && eds->update_calls == eds->update_period) {
    

    //use estimated variance to take a step
    real step_size = 0;
    real temp;

    for(i = 0; i < eds->cv_number; i++) {
      //calulcate step size
      temp = 2. * (eds->means[i] / eds->centers[i] - 1) * eds->ssd[i] / 
	(eds->update_calls - 1);
      step_size = temp / (eds->simtemp * boltz);

      //reset means/vars
      eds->means[i] = 0;
      eds->ssd[i] = 0;
      
      //multidimesional stochastic step
      if(eds->cv_number == 1 || rando(&eds->seed) < 1. / eds->cv_number) {
	eds->coupling_accum[i] += step_size * step_size;
	eds->current_coupling[i] = eds->set_coupling[i];
	eds->set_coupling[i] += eds->max_coupling_range[i] / 
	  sqrt(eds->coupling_accum[i]) * step_size;
	eds->coupling_rate[i] = (eds->set_coupling[i] - eds->current_coupling[i]) / eds->update_period;
	eds->coupling_rate[i] = copysign(fmin(fabs(eds->coupling_rate[i]),
					     eds->max_coupling_rate[i]), 
					eds->coupling_rate[i]);
	
      } else {
	//we chose not to change the bias
	eds->coupling_rate[i] = 0;
      }      
    } // closing colvar loop
    
    eds->update_calls = 0;
    eds->b_equilibration = true; //back to equilibration now
  } //close if update if

  return bias_energy;
}
