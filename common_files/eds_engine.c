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
		     const char* filename,
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

  eds->output_filename = (char*) malloc(sizeof(char) * (strlen(filename) + 1));
  strcpy(eds->output_filename, filename);
  eds->output_file = NULL;

  //divide it by 2 so we spend half equilibrating and half collecting statistics.
  eds->update_period = update_period / 2;
  eds->update_calls = 0;
  eds->b_equilibration = 1;
  eds->b_hard_coupling_range = b_hard_coupling_range;

  int i;
  for(i = 0; i < cv_number; i++) {
    eds->max_coupling_range[i] = 1;
  }
 
}

void PREFIX eds_read(char **word, int nw, t_plumed_input *input, FILE *fplog) {
  //EXAMPLE for adaptive
  //EDS STRIDE 500 SIMTEMP 300 SEED 4143 FILENAME FOO CV LIST 1 2 4
  //EDS CV CENTERS 0.5 2.5 2.3
  //EDS CV RANGES 5 5 5
  //[optional] EDS CV CONSTANTS 2403.1 4003.31 499.1
  //EXAMPLE for fixed
  //EDS SIMTEMP 300 CV LIST 1 2 4
  //EDS CV CENTERS 0.5 2.5 2.3
  //EDS CV RANGES 34 23 4      
  //EXAMPLE for restart
  //EDS STRIDE 500 SIMTEMP 300 SEED 4143 FILENAME FOO RESTART BAR CV LIST 1 2 4
  //EDS CV CENTERS 0.5 2.5 2.3
  //EDS CV RANGES 5 5 5

  


  int i, icv, iw, iat, j;
  real uno;
  char restart_filename[200];
  int update_period = 0;
  int eds_seed = 0;
  char filename[200];
  int restart = 0;
  int* cv_map = NULL;

  if(!logical.eds) 
    fprintf(fplog, "Enabling experiment directed simulation\n");
  logical.eds = 1;
  iw = 1;
  
  if(!strcmp(word[iw], "CV")) {
    iw++;
    if(!strcmp(word[iw++], "CENTERS")) {
      if(eds.cv_number == 0) {
	plumed_error("Must define CVs first for EDS with [CV LIST 1 2 3]");
      }
      for(icv = 0;iw < nw; iw++) {
	sscanf(word[iw], "%lf", &uno);
	eds.centers[icv] = uno;
	fprintf(fplog, "EDS: Will center CV %d at %lf\n", icv+1,
		eds.centers[icv]);
	icv++;
      }
    }
    else if(!strcmp(word[iw - 1], "RANGES")) {
      if(eds.cv_number == 0) {
	plumed_error("Must define CVs first for EDS with [CV LIST 1 2 3]");
      }
      for(icv = 0;iw < nw; iw++) {
	sscanf(word[iw], "%lf", &uno);
	eds.max_coupling_range[icv] = uno;
	eds.max_coupling_rate[icv] = eds.max_coupling_range[icv] / (10 * eds.update_period);
	fprintf(fplog, 
		"EDS: Will cap range of CV %d at %lf\n", 
		icv+1,
		eds.max_coupling_range[icv]);
	icv++;
      }
    } else if(!strcmp(word[iw - 1], "CONSTANTS")) {
      if(eds.cv_number == 0) {
	plumed_error("Must define CVs first for EDS with [CV LIST 1 2 3]");
      }
      for(icv = 0;iw < nw; iw++) {
	sscanf(word[iw], "%lf", &uno);
	eds.current_coupling[icv] = uno;
	eds.set_coupling[icv] = uno;
	fprintf(fplog, 
		"EDS: Starting CV %d at %lf\n", 
		icv+1,
		eds.current_coupling[icv]);
	icv++;
      }
    } else {
      plumed_error("Syntax is EDS CV RANGES.... or EDS CV CENTERS....\n");
    }
  } else {
    
    if(eds.cv_number != 0) {
      plumed_error("Syntax is EDS CV RANGES.... or EDS CV CENTERS....or EDS CV CONSTANTS\n");
    }

    int last_iw;
    uno = -1;
    while(iw < nw) {
      last_iw = iw;
      if(!strcmp(word[iw], "STRIDE"))
	if(!sscanf(word[++iw], "%d", &update_period))
	  plumed_error("Could not read stride\n");	  
	else
	  iw++;
      
      if(!strcmp(word[iw], "SIMTEMP")) {
	if(!sscanf(word[++iw], "%lf", &uno)){
	  plumed_error("Must specify SIMTEMP in EDS [EDS STRIDE 500 SIMTEMP 300 SEED 4313 CV LIST 1 3]\n");
	}
	iw++;
      } else if(uno == -1) {
	plumed_error("Must specify SIMTEMP in EDS as first argument [EDS STRIDE 500 SIMTEMP 300 SEED 431 CV LIST 1 3]\n");
      }
      
      if(!strcmp(word[iw], "SEED")) 
	if(!sscanf(word[++iw], "%d", &eds_seed))
	  plumed_error("Must use integer SEED\n");
	else
	  iw++;
      
      if(!strcmp(word[iw], "FILENAME")) {
	if(!sscanf(word[++iw], "%s", filename)){
	  plumed_error("Filename invalid\n");
	} else {
	  iw++;
	}
      } else {
	strcpy(filename, "EDS_OUT");
      }
      
      if(!strcmp(word[iw], "RESTART")) {
	if(!sscanf(word[++iw], "%s", restart_filename)){
	  plumed_error("Restart Filename invalid\n");
	} else {
	  iw++;
	  restart = 1;
	}
      }
      if(!strcmp(word[iw], "CV") && !strcmp(word[++iw], "LIST")) {
	iw++;
	cv_map = (int*) malloc(sizeof(int) * nconst_max);
	for(i = 0;iw < nw; iw++) {
	  sscanf(word[iw], "%d", &icv);
	  cv_map[i] = icv - 1;
	  fprintf(fplog, 
		  "EDS: Will use CV %d \n", 
		  icv);
	  i++;
	}
      }
      if(last_iw == iw) {
	fprintf(fplog,"WARNING: Ignoring invalid option %s\n", word[iw]);
	iw++;
      }
    }
    if(cv_map == NULL)
      plumed_error("Must specify CV List in EDS as last argument [EDS 500 CV LIST 1 3]\n");
    cv_map = (int *) realloc(cv_map, sizeof(int) * i);
    eds_init(i, update_period, uno, eds_seed, 0, cv_map, (const char*) filename, &eds);   
    if(restart)
      eds_read_restart(restart_filename, fplog, &eds);
  }
}

void PREFIX eds_read_restart(char* filename, FILE* fplog, t_eds* eds) {

  int i;

  //Just read it over and over again, so that only the last line is used
  FILE* restart = fopen(filename, "r");
  if(restart == NULL) {
    //just skip with warning, since it might be included for future restarts
    fprintf(fplog, "WARNING: Could not find restart file in EDS, skipping...\n");
  }
  int success = 1;
  long long int temp;

  fprintf(fplog, "READING IN EDS RESTART...");
  while(!feof(restart)) {
    //read and skip step
    success &= fscanf(restart, "%lld ", &temp);
    for(i = 0; i < eds->cv_number; i++) {
      success &= fscanf(restart, "%lf", &(eds->current_coupling[i]));
    }
    for(i = 0; i < eds->cv_number; i++) {
      success &= fscanf(restart, "%lf ", &(eds->coupling_accum[i]));    
    }
    if(!success)
      plumed_error("Found incomplete line in EDS restart file");
  }

  fprintf(fplog, "DONE\n");

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
  free(eds->output_filename);
  
  fclose(eds->output_file);
  

}

/*
 * Run update of calculation and apply forces
 */
real PREFIX eds_engine(real* ss0, real* force, 
		       t_eds* eds, real boltz) {

  real bias_energy = 0.0;
  int i;  


  if(eds->update_calls == 0 && eds->update_period != 0) {
    for(i = 0; i < eds->cv_number; i++) {
      eds->max_coupling_rate[i] = eds->max_coupling_range[i] / (10 * eds->update_period);
    }
  }

  eds->update_calls++;

  int b_finished_equil_flag = 1;
  real delta;

  //zero forces
  for(i = 0; i < eds->cv_number; i++)
    force[eds->cv_map[i]] = 0;

  //apply forces for this setp and calculate energies
  for(i = 0; i < eds->cv_number; i++) {
    force[eds->cv_map[i]] -= eds->current_coupling[i] / eds->centers[i];
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

void PREFIX dump_array(real* array, int length, FILE* file, const char* name) {

  int i;
  fprintf(file, "%s: ", name);
  for(i = 0; i < length; i++) {
    fprintf(file, "%0.3f ", array[i]); 
  }
  fprintf(file, "\n");

}

void PREFIX eds_dump(t_eds* eds) {
  
  dump_array(eds->centers, eds->cv_number, eds->output_file, "centers");
  dump_array(eds->means, eds->cv_number, eds->output_file, "means");
  dump_array(eds->ssd, eds->cv_number, eds->output_file, "ssd");
  dump_array(eds->max_coupling_range, eds->cv_number, eds->output_file, "max_coupling_range");
  dump_array(eds->max_coupling_rate, eds->cv_number, eds->output_file, "max_coupling_rate");
  dump_array(eds->set_coupling, eds->cv_number, eds->output_file, "set_coupling");
  dump_array(eds->current_coupling, eds->cv_number, eds->output_file, "current_coupling");
  dump_array(eds->coupling_rate, eds->cv_number, eds->output_file, "coupling_rate");
  dump_array(eds->coupling_accum, eds->cv_number, eds->output_file, "coupling_accum");
  
  int i;
  fprintf(eds->output_file, "%s: ", "cv_map");
  for(i = 0; i < eds->cv_number; i++) {
    fprintf(eds->output_file, "%5d ", eds->cv_map[i]); 
  }
  fprintf(eds->output_file, "\n");

  fprintf(eds->output_file, "simtemp: %f\n", eds->simtemp);
  fprintf(eds->output_file, "cv_number: %d\n", eds->cv_number);
  fprintf(eds->output_file, "update_period: %d\n", eds->update_period);
  fprintf(eds->output_file, "update_calls: %d\n", eds->update_calls);

}


void PREFIX eds_write(t_eds* eds, long long int step) {

  if(eds->update_period > 0 && eds->update_calls % eds->update_period == 0) {

    if(eds->output_file == NULL)
      eds->output_file = fopen(eds->output_filename, "w");
    
    int i;
    
    fprintf(eds->output_file, "%12ld ", step);
    
#ifndef DUMP_EDS
    
    for(i = 0; i < eds->cv_number; i++)
      fprintf(eds->output_file, "%0.5f ", eds->current_coupling[i]);
    for(i = 0; i < eds->cv_number; i++)
      fprintf(eds->output_file, "%0.5f ", eds->coupling_accum[i]);
    
    //flush file
    fprintf(eds->output_file, "\n");
    fflush(eds->output_file);
#else
    fprintf(eds->output_file, "\n");
    eds_dump(eds);
    fprintf(eds->output_file, "-------------------------");
#endif//DUMP_EDS
  }
  
}
