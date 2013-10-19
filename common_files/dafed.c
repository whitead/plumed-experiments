/*
* dafed.c
*
*  This file is an addition to plumed for d-AFED
*      Author: Michel Cuendet
*
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
#include <math.h>



// =====================================================================================
void PREFIX initialize_ggmt(
		struct ggmt_s *ggmt,				// Generalized Gaussian Moments Thermostat structure
		double T0,						// reference temperature
		double tau,					// time constant
		double  delta_t				// external MD time step
	) {

	/* We assume that temperature and tau and n_respa are given */
	ggmt->kT0	= mtd_data.boltz * T0;
	ggmt->Q1	= ggmt->kT0 * tau * tau;
	ggmt->Q2	= 8.0/3.0 * pow(ggmt->kT0,3) *tau *tau ;

	ggmt->dt[0]  = delta_t /((double)ggmt->n_respa_ggmt);
   /* Weights for the Suzuki-Yoshida decomposition : */
   /* These correspond to the delta t_j from Liu et al. (2000) */
	ggmt->dt[1]  = (1.0 / (2.0 - pow(2.0,(1.0/3.0)))) * ggmt->dt[0]  ;
	ggmt->dt[2]  = ggmt->dt[0] -  2.0*ggmt->dt[1];
	ggmt->dt[3]  = ggmt->dt[1];

	ggmt->eta1	= 1.0;
	ggmt->eta2	= 1.0;
	ggmt->v1	= -sqrt( ggmt->kT0 / ggmt->Q1 );
	ggmt->v2	= -sqrt( ggmt->kT0 / ggmt->Q2 );


}

// =====================================================================================
void PREFIX initialize_dafed(struct dafed_s *daf, real dt)
	{ 												// dt is taken as real from mtd_data

	  // This is called after reading the input.
	  // daf.s is initialized later to the value of colvar.ss0[i_c]
	  daf->vs = -sqrt(mtd_data.boltz * daf->temperature / daf->mass );
	  daf->dt_respa = dt;  // In the current implementation, RESPA is done in the MD loop.
	  daf->f = 0.0;
	  daf->do_initialize_s = 1;
	  daf->do_skip_integration = 0;
	  daf->thermo_work = 0.0;
	  daf->dafed_work = 0.0;
	  daf->do_jacobian_force = 0;
	  daf->jacobian_force = 0.0;
	  daf->dWold = 0.0;

	  initialize_ggmt(&(daf->ggmt), daf->temperature, daf->tauthermo, daf->dt_respa);

}

// =====================================================================================
void PREFIX get_kt(struct dafed_s *daf, double *d) {
  *d = pow(daf->vs,2) * daf->mass;
}

// =====================================================================================
void PREFIX integrate_ggmt(struct dafed_s *daf) {

		// This corresponds to the integration over HALF an external time step
		// This routine has to be called twice per step

       struct ggmt_s *th = &(daf->ggmt);
	   int iii, jjj;
	   double aa=0.0, bb=0.0;
	   double kt;
	   double G1,G2;
	   double dt2, dt4, dt8;

	   get_kt(daf, &kt);				   // kt is equivalent to (p^2 / m)

	   // Loop over the RESPA integration steps
	   for (iii = 1; iii <= th->n_respa_ggmt; iii++) {
		 // Loop over the three therms of the Suzuki-Yoshida decomposition :
	     for (jjj = 1; jjj <= 3; jjj++) {

	       dt2 = 0.5*th->dt[jjj];   // half time step
	       dt4 = 0.25*th->dt[jjj];	// quarter time step
	       dt8 = 0.125*th->dt[jjj]; // 1/8 time step

	       G1   = (kt - th->kT0)/(th->Q1);
	       G2   = (pow(kt,2)/(3.0) - pow(th->kT0,2))/(th->Q2);
	       th->v1 += dt4 * G1;
	       th->v2 += dt4 * G2;

	       aa         = exp(-dt8 * (th->v1 + th->kT0*th->v2));
	       daf->vs   *= aa;
	       get_kt(daf, &kt);

	       bb         = kt*th->v2/(3.0);
	       daf->vs   *= sqrt(1.0/(1.0 + dt2 * bb));

	       daf->vs   *= aa;
	       get_kt(daf, &kt);

	       th->eta1 += dt2 * th->v1;
	       th->eta2 += dt2 * th->v2*(th->kT0 + kt);

	       aa         = exp(-dt8 * (th->v1 + th->kT0*th->v2));
	       daf->vs   *= aa;
	       get_kt(daf, &kt);

	       bb         = kt*th->v2/(3.0);
	       daf->vs   *= sqrt(1.0/(1.0 + dt2 * bb));

	       daf->vs   *= aa;
	       get_kt(daf, &kt);

	       G1   = (kt - th->kT0)/(th->Q1);
	       G2   = (pow(kt,2)/(3.0) - pow(th->kT0,2))/(th->Q2);
	       th->v1 += dt4 * G1;
	       th->v2 += dt4 * G2;

	     }


	   }
}

// =====================================================================================
real PREFIX energy_ggmt(struct ggmt_s ggmt) {
	real ener;

	/* In the cas of DAFED, number of particles =1, dimention = 1 always */
	ener = 0.5* (  ggmt.Q1*pow(ggmt.v1,2)  +
				   ggmt.Q2*pow(ggmt.v2,2) )  +
				   ggmt.kT0*(ggmt.eta1 + ggmt.eta2);
	return ener;
}

// =====================================================================================

void PREFIX print_dafed(struct dafed_s *daf,FILE *cv_file, int i_c)
  {
		double potential;
		double energy;
		double temperature;
		double transferred_work;
		double delta;

		delta =  daf->colvar - daf->s;
		if (daf->do_periodicity) {
				if (delta > 0.5 * daf->periodicity_gap ) delta = delta - daf->periodicity_gap;
				if (delta < - 0.5 * daf->periodicity_gap ) delta = delta + daf->periodicity_gap;
		}

		potential =  0.5*daf->kappa * pow(delta ,2);
		energy = 0.5*daf->mass * pow(daf->vs,2) +
						  potential +
						  energy_ggmt(daf->ggmt);

		temperature = pow(daf->vs,2) * daf->mass / mtd_data.boltz;

		// What we want is not the accumulated work which pertains to the system H_phys + V
		// But we want the work W transferred to the physical system, which pertains to H_phys
		// Where H_phys is the energy of the physical system + thermostat
		// and V is the potential 0.5*kappa*(q-s)^2
		// For a discussion of this (in the context of steered MD), see
		// Schurr and Fujimoto, J. Phys. Chem. B 107. 14007 (2003)
		// This way, both H-W and H_dafed+W should be conserved
		transferred_work = daf->dafed_work - potential;

  		 // fprintf(cv_file,"      d-AFED %d ",i_c+1);
  		 fprintf(cv_file," %10.5e",daf->s);
  		 fprintf(cv_file," %10.5f",temperature);
  		 fprintf(cv_file," %10.5e",energy);
  		 fprintf(cv_file," %13.8e",transferred_work);
  		 // fprintf(cv_file," %15f",daf->thermo_work);

  }

// =====================================================================================
void PREFIX read_dafed_state()
{
		int i_c;
		FILE *file;
		double time;


		if (!logical.do_dafed)
			plumed_error("CANNOT READ d-AFED STATE IF NO d-ADFED CV IS DEFINED");

		file = fopen(dafed_control.in_file, "r");
		if(!file) {										// not exist
		           char buf[1024];
		           sprintf(buf,"Cannot read d-AFED state file %s for reading !",dafed_control.in_file);
		           plumed_error(buf);
		   }

		fscanf(file, "%lf\n",&time);

		for(i_c=0;i_c<colvar.nconst;i_c++){
			  if (logical.dafed[i_c]) {
				    fscanf(file, "%lf %lf %lf %lf %lf %lf\n",
				    		&(dafed[i_c].s),
				    		&(dafed[i_c].vs),
				    		&(dafed[i_c].ggmt.eta1),
				    		&(dafed[i_c].ggmt.v1),
				    		&(dafed[i_c].ggmt.eta2),
				    		&(dafed[i_c].ggmt.v2));

				    fprintf(mtd_data.fplog, "|- DAFED CV %d : READ RESTART VALUES \n",i_c+1);
				    fprintf(mtd_data.fplog, "|- \t s \t %.15e\n",dafed[i_c].s);
				    fprintf(mtd_data.fplog, "|- \t vs \t %.15e\n",dafed[i_c].vs);
				    fprintf(mtd_data.fplog, "|- \t eta1 \t %.15e\n",dafed[i_c].ggmt.eta1);
				    fprintf(mtd_data.fplog, "|- \t v1 \t %.15e\n",dafed[i_c].ggmt.v1);
				    fprintf(mtd_data.fplog, "|- \t eta2 \t %.15e\n",dafed[i_c].ggmt.eta2);
				    fprintf(mtd_data.fplog, "|- \t v2 \t %.15e\n",dafed[i_c].ggmt.v2);
			  }
		}
		fclose(file);
}

// =====================================================================================
void PREFIX write_dafed_state()
{
		int i_c;
		char stringa[500];
		FILE *file;

		sprintf(stringa, "%s.old", dafed_control.out_file);
		if(mtd_data.ionode) rename(dafed_control.out_file, stringa);

		file = fopen((mtd_data.ionode?dafed_control.out_file:"/dev/null"), "w");
		// Here we use the same construct as for the COLVAR file, where all nodes write to /dev/null except the ionode.
		if(!file) {										// not exist
		           char buf[1024];
		           sprintf(buf,"Cannot open d-AFED state file %s for writing !",dafed_control.out_file);
		           plumed_error(buf);
		   };

		fprintf(file, "%10.3f\n", mtd_data.time);
		for(i_c=0;i_c<colvar.nconst;i_c++){
			  if (logical.dafed[i_c]) {
				    fprintf(file, "%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e\n",
				    		dafed[i_c].s,
				    		dafed[i_c].vs,
				    		dafed[i_c].ggmt.eta1,
				    		dafed[i_c].ggmt.v1,
				    		dafed[i_c].ggmt.eta2,
				    		dafed[i_c].ggmt.v2);
			  }
		}
		fclose(file);
}


// =====================================================================================
void PREFIX dafed_engine(real *this_colvar)
{											  // this_colvar is taken as real from colvar.ss0
	int i,i_c;
	double kt;
	double dW;
	double delta;
	struct dafed_s *daf;
	int do_write=0;

	for(i_c=0;i_c<colvar.nconst;i_c++){
	  if (logical.dafed[i_c]) {

		    daf = &dafed[i_c];

		    if (daf->do_initialize_s) {
				daf->s = (double)this_colvar[i_c];
				daf->do_initialize_s =0;
			}

			daf->colvar = (double)this_colvar[i_c];

			if (daf->do_skip_integration) {
				// This is if we have just read a restart file
				daf->f = - daf->kappa * ( daf->colvar  - daf->s);
				daf->do_skip_integration = 0;
			}else{

			  //for(i=1;i<=daf->n_respa;i++) {
			  // In the current implementation, RESPA is done in the MD loop

			  // integrate from 0 to dt/2
				integrate_ggmt(daf);
				// We want the force FROM the colvar ON s -> change sign
				daf->vs +=  0.5 * daf->dt_respa * (-daf->f + daf->jacobian_force) / daf->mass ;
				daf->s +=  daf->dt_respa * daf->vs;
				if (daf->do_periodicity) {
					if (daf->s < daf->periodicity_low) daf->s = daf->s + daf->periodicity_gap;
					if (daf->s > daf->periodicity_high) daf->s = daf->s - daf->periodicity_gap;
				}

			  // Integrate the thermostat work : (for now only a simple rectangular scheme)
				// In principle here we have kt and v1, v2 at dt/2
				get_kt(daf, &kt);
				daf->thermo_work -= 0.5*daf->dt_respa  *
							kt*( daf->ggmt.v1 +  daf->ggmt.v2*(daf->ggmt.kT0 + kt/3.0) );

			  // calculate the force.
				// This is the DAFED force at full time step t
				// (i.e. after the last RESPA iteration is completed)
				// This is the force FROM s ON the colvar -> minus sign
				delta = daf->colvar  - daf->s ;
				if (daf->do_periodicity) {
					if (delta > (0.5 * daf->periodicity_gap) ) delta = delta - daf->periodicity_gap;
					if (delta < (- 0.5 * daf->periodicity_gap) ) delta = delta + daf->periodicity_gap;
				}
				daf->f = - daf->kappa * delta;

				if ( daf->do_jacobian_force ) {
					daf->jacobian_force = - 2.0 * 300.0 * 0.00831451 / daf->s ;
				}

			  // integrate from dt/2 to dt
				// We want the force FROM the colvar ON s -> change sign
				daf->vs +=  0.5 * daf->dt_respa * (-daf->f + daf->jacobian_force) / daf->mass ;
				integrate_ggmt(daf);

			  // Integrate the thermostat work : (for now only a simple rectangular scheme)
				// In principle here we have kt and v1, v2 at full dt
				get_kt(daf, &kt);
				daf->thermo_work -= 0.5*daf->dt_respa  *
								kt*( daf->ggmt.v1 +  daf->ggmt.v2*(daf->ggmt.kT0 + kt/3.0) );

			 // Integrate DAFED work : (triangle scheme)
				dW = daf->dt_respa * daf->f * daf->vs;
				daf->dafed_work += 0.5 * ( dW + daf->dWold);
				daf->dWold = dW;
			  //}
			} // if skip
	  } // if
	} // for

	if (dafed_control.write_freq == -1)
		do_write = (logical.not_same_step)&&(!firstTime)&&(dafed_control.do_cpt);
	else if (dafed_control.write_freq > 0)
		do_write = (logical.not_same_step)&&(!firstTime)&&(!(colvar.it%dafed_control.write_freq));

	if (do_write) write_dafed_state();

}


