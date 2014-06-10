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

//ADW>

int independent_cv_cache_natoms;
int independent_cv_cache_natoms2;
int independent_cv_cache_atom;
int independent_cv_cache_atom2;

/*
 * Makes a CV only have one atom. Returns 1 if the given atom_index is
 * valid. It returns 0 if it should be called again and -1 if finished
 */
int PREFIX independent_stash_cv(int i_c, int atom_index){

  switch(colvar.type_s[i_c]){ 
  case 32:
    //restraint_position 
    if(atom_index < colvar.natoms[i_c]) {
      //replace atom number and position with cached version
      independent_cv_cache_natoms = colvar.natoms[i_c];
      colvar.natoms[i_c] = 1;
      independent_cv_cache_atom = colvar.cvatoms[i_c][0];
      colvar.cvatoms[i_c][0] = colvar.cvatoms[i_c][atom_index];
      return 1;
    }
    break;
  case 1:
    //restraint_dist
    // colvar.list is length of list 1. colvar.natoms is length of list 1 + 2

    //check if we're done
    if(atom_index < (colvar.natoms[i_c] - colvar.list[i_c][0]) * colvar.list[i_c][0]) {

      int index1 = atom_index % colvar.list[i_c][0];
      int index2 = colvar.natoms[i_c] - colvar.list[i_c][0]
	+ atom_index / colvar.list[i_c][0];

      //Check if the pair is the same particle
      if(colvar.cvatoms[i_c][index1] == colvar.cvatoms[i_c][index2])
	return 0;
    
      //this is a pair-wise list that needs to be overwritten so one pair
      //is processed at a time        
      //first pair
      independent_cv_cache_atom = colvar.cvatoms[i_c][0];
      colvar.cvatoms[i_c][0] = colvar.cvatoms[i_c][index1];
      //second pair
      independent_cv_cache_atom2 = colvar.cvatoms[i_c][1];
      colvar.cvatoms[i_c][1] = colvar.cvatoms[i_c][index2];

      
      //store the length of the two lits
      independent_cv_cache_natoms = colvar.natoms[i_c];
      independent_cv_cache_natoms2 = colvar.list[i_c][0];
      
      //make lengths 1
      colvar.list[i_c][0] = 1;
      colvar.natoms[i_c] = 2;

      return 1;
    }
    break;
  }
  return -1;
}

/*
 * Undoes the changes from independent_stash_cv. Returns 1 if
 * stash was done, and thus forces should be added.
 */
int PREFIX independent_pop_cv(int i_c, int atom_index){
   //swap values with the cache

  switch(colvar.type_s[i_c]) { 
  case 32:
    //restraint_position 
    if(colvar.natoms[i_c] == 1) {
      colvar.natoms[i_c] = independent_cv_cache_natoms;
      colvar.cvatoms[i_c][atom_index] = independent_cv_cache_atom;
      return 1;
    }
    break;
  case 1:
    //restraint_position
    if(colvar.list[i_c][0] == 1) {
      //first pair
      colvar.cvatoms[i_c][0] = independent_cv_cache_atom;
      //second pair
      colvar.cvatoms[i_c][1] = independent_cv_cache_atom2;
      //store the length of the two lits
      colvar.natoms[i_c] = independent_cv_cache_natoms;
      colvar.list[i_c][0] = independent_cv_cache_natoms2;
      return 1;
    }
    break;
  }
  return 0;
}

//<ADW

void PREFIX restraint(struct mtd_data_s *mtd_data)
{
  int i_c, ind_i_c, ncv, nth, ntrh, ntp, rpxm, ntwg;                    // indexes for cycles and time

  hills.Vhills = cvw.Vwall = Vext = Vrecon = Vconstr = Vsteerplan = 0.;  // Hills,  Walls , Constraint, external potential and steerplan energy initialization

  colvar.it=mtd_data->istep;                                   // update microdynamics step
  if(colvar.it==mtd_data->istep_old){
    logical.not_same_step=0;
  } else {
    logical.not_same_step=1;
    mtd_data->istep_old=colvar.it;
  }

  ncv  = colvar.nconst;                                                                           	                  // number of CVs
  ntp  = (logical.not_same_step)&&(logical.print)&&(!(colvar.it%colvar.nt_print));			                  // have I got to print COLVAR?
  nth  = ( (logical.not_same_step)&&(logical.do_hills)&&(!(colvar.it%hills.nt_hills))&&(!firstTime) )
         || (hills.max_height>0);              									          // have I got to add HILL?
  ntrh = (logical.not_same_step)&&(logical.do_hills)&&(!(colvar.it%hills.nr_hills))&&(logical.do_walkers)&&(!firstTime) ; // period in steps to read HILLS
  ntwg = (logical.write_grid)&&(logical.not_same_step)&&(!(colvar.it%bias_grid.w_stride))&&(!firstTime);                       // write GRID on file

#ifdef PLUMED_GROMACS
  set_pbc(&mtd_data->metapbc, mtd_data->ePBC, mtd_data->cell);                                              // set pbc
#endif



  if(logical.rpxm) rpxm = !(colvar.it%mtd_data->repl_ex_nst) && (logical.not_same_step);                    // have I got a replica exchange trial
  else rpxm = 0;
  if(logical.debug){						// debug is a GROMACS global which identifies
    test_derivatives(mtd_data);		                       // the use of -debug options in mdrun
    #ifdef RECONMETAD
     #ifndef DRIVER
       if( reconOn==1 ){test_recon_derivatives(mtd_data);}    // Reconnaissance metadynamics version of test derivatives
     #endif
    #endif  
    EXIT();			        // allow for only one step of dynamics
  }

// eventually align atoms
  if(colvar.align_atoms){
     int i,iatom1,iatom2;
     real distance[3],dummy;
     for(i=0;i<colvar.align_atoms-1;i++){
       iatom1=colvar.align_list[i];
       iatom2=colvar.align_list[i+1];
       minimal_image(mtd_data->pos[iatom2],mtd_data->pos[iatom1],&dummy,distance);
       mtd_data->pos[iatom2][0]=mtd_data->pos[iatom1][0]+distance[0];
       mtd_data->pos[iatom2][1]=mtd_data->pos[iatom1][1]+distance[1];
       mtd_data->pos[iatom2][2]=mtd_data->pos[iatom1][2]+distance[2];
     }
   };

// eventually dump atoms
#ifdef PLUMED_GROMACS
  if(mtd_data->dump_atoms && (logical.not_same_step)&&(!(colvar.it%mtd_data->dump_stride)) ){
    rvec* atoms;
    int i,iat;
    if(!mtd_data->dump_file) mtd_data->dump_file = open_xtc(mtd_data->dump_filen,"w");
    snew(atoms,mtd_data->dump_atoms);
    for(i=0;i<mtd_data->dump_atoms;i++){
      iat=mtd_data->dump_list[i];
      atoms[i][0]=mtd_data->pos[iat][0];
      atoms[i][1]=mtd_data->pos[iat][1];
      atoms[i][2]=mtd_data->pos[iat][2];
    };
    write_xtc(mtd_data->dump_file,mtd_data->dump_atoms,colvar.it,mtd_data->time,mtd_data->cell,atoms,1000.0);
    sfree(atoms);
  };
#endif

  //ADW>
  //repeat the entire algorithm for each independent CV   
  int stash_result = 1;
  int remaining_ind = 1;
  //Zero forces, moved here otherwise we overwrite the forces for each independent CV loop iteration
  zero_forces(mtd_data);

  for(ind_i_c = 0; ind_i_c < remaining_ind; ind_i_c++) {

    //try stash and see if more are necessary
    for(i_c=0;i_c<ncv;i_c++) {
      if(colvar.b_treat_independent[i_c]) {
	stash_result = independent_stash_cv(i_c, ind_i_c);//try to stash the other atoms and use only ind_i_c
	if(stash_result != 1)
	  break;
      }
    }

    //pop stash if we didn't succeed
    if(stash_result != 1)  {
      for(i_c=0;i_c<ncv;i_c++)
	if(colvar.b_treat_independent[i_c])
	  independent_pop_cv(i_c, ind_i_c);

      if(stash_result == 0) {
	remaining_ind++;
	continue; //try again
      }
      break; //we're done, stash_result == -1
    } else
      remaining_ind++; //make sure we try another iteration

    //<ADW

    // this cycle is intended to calculate CVs values and derivatives

    for(i_c=0;i_c<ncv;i_c++){
      colvar.ff_hills[i_c] = 0.;                 	// initialization hills forces
      cvw.fwall[i_c]       = 0.;  		// initialization walls forces
      fext[i_c]            = 0.;                  // initialization external forces
      
      //    if((!logical.always[i_c])&&(logical.rpxm)&&(!rpxm)&&(!ntp)) continue;
      if((!logical.always[i_c])&&(!ntp)) {   //if a variable is used only for printing purposes and it isn't the time to print then
	if(logical.rpxm) {
	  if((!rpxm)&&(!firstTime)) continue;  //if we are doing bias-exchange and it isn't the time for an exchange we could avoid to calculate the variable
	} else continue; //if we are not doing bias-exchange we could avoid to calculate the variable
	// but if we are doing bias-exchange and it is time to try an axchange we must calculate all the variables
      }
      
      //fprintf(mtd_data->fplog,"CALCULATING CV %d\n",i_c); 
      //fflush(mtd_data->fplog); 
      
      switch(colvar.type_s[i_c]){
	// geometric CVs
      case 1: dist_restraint(i_c, mtd_data); break;			// DISTANCE
      case 2: mindist_restraint(i_c, mtd_data); break;               	// MINDIST
      case 3: coord_restraint(i_c, mtd_data); break;	              	// COORD
      case 4: angle_restraint(i_c, mtd_data); break;	                // ANGLE
      case 5: torsion_restraint(i_c, mtd_data); break;                  // TORSION
      case 6: alfabeta_restraint(i_c, mtd_data); break;                	// ALPHA-BETA
	// interaction CVs
      case 7: hbonds_restraint(i_c, mtd_data); break;                   // HBONDS
      case 8: dipole_restraint(i_c, mtd_data); break;       		// DIPOLE
	// conformations CVs
      case 11: radgyr_restraint(i_c, mtd_data); break;   	       	// RGYR
      case 16: dihcor_restraint(i_c, mtd_data); break;                 	// DIHEDRAL-COR
	// water CVs
      case 20: waterbridge_restraint(i_c, mtd_data); break;            	// WATERBRIDGE
      case 21: density_restraint(i_c, mtd_data); break;                 // DENSITY
      case 22: densityswitch_restraint(i_c, mtd_data); break;           // DENSITYSWITCH
	// trajectory CVs
      case 30: spath_restraint(i_c, mtd_data); break;                   // S_MAPPATH
      case 31: zpath_restraint(i_c, mtd_data); break;                   // Z_MAPPATH
      case 32: position_restraint(i_c, mtd_data); break;                // ATOM POSITION
      case 33: elstpot_restraint(i_c, mtd_data); break;                 // ELSTPOT POSITION
      case 34: puckering_restraint(i_c, mtd_data); break;               // PUCKERING
      case 35: energy_restraint(i_c, mtd_data); break;                  // ENERGY
      case 36: helix_restraint(i_c, mtd_data); break;                   // HELIX
      case 37: alpharmsd_restraint(i_c, mtd_data); break;               // ALPHARMSD
      case 38: antibetarmsd_restraint(i_c, mtd_data); break;            // ANTIBETARMSD
      case 39: parabetarmsd_restraint(i_c, mtd_data); break;            // PARABETARMSD
	//case 40: camshift_restraint(i_c, mtd_data); break;              // CAMSHIFT ENERGY
	//case 41: camshiftens_restraint(i_c, mtd_data); break;           // ENS CAMSHIFT ENERGY
      case 42: pca_restraint(i_c, mtd_data); break;      		// PCA PROJECTION
      case 45: cmap_restraint(i_c, mtd_data); break;                    // CMAP 
#if defined CVS  
      case 46: bespoke_restraint(i_c, mtd_data); break;                 // BESPOKE COLLECTIVE COORDINATES
#endif
      case 47: rdf_restraint(i_c, 1, mtd_data); break;                  // DISCRETIZED RDF
      case 49: histogram_restraint(i_c, mtd_data); break;               // HISTOGRAM CV
      case 50: poly_restraint(i_c, mtd_data); break;			// POLYNOMIAL COMBINATION
      case 51: func_restraint(i_c, mtd_data); break;			// GENERAL FUNC OF CVS 
      case 52: adf_restraint(i_c, 1, mtd_data); break;                  // DISCRITIZED ANGLE DISTRIBUTION FUNCTION 
      case 53: msd_restraint(i_c,  mtd_data); break;    	
      case 55: sprint_restraint(i_c, mtd_data); break;                  // SPRINT TOPOLOGICAL CV
      }
      
#ifdef PATHREF_FINDIFF
      fprintf(mtd_data->fplog,"|---CALLING THE TEST \n");
      //finite difference tests over reference frame derivatives
      if(colvar.type_s[i_c]==30){pathref_findiff(i_c,mtd_data);}
      if(colvar.type_s[i_c]==31){pathref_findiff(i_c,mtd_data);}
      fprintf(mtd_data->fplog,"|---END OF CALL \n");
      EXIT();
#endif
     
      //ADW>
      //handle independently treated CVs
      if(colvar.b_treat_independent[i_c])
	independent_pop_cv(i_c, ind_i_c);
      //<ADW
 
    }
    
    mtd_data->time=colvar.it*(mtd_data->dt)+mtd_data->time_offset;
    
    if(logical.commit) commit_analysis();     // committors analysis
    
    if(couplingmatrix.is_on) calc_couplingmatrix(colvar.it);	
    
    // this is the really dynamics code in which we calculate hills forces and then forces on CVs.
    if(logical.do_hills){
      hills_force();						// compute hills force and energy
      if(logical.widthadapt) hills_adapt();			// to adapt gaussian width
      if(nth) {
	//if we are stochastically sampling, do it now
	if(colvar.stoch_sample == 1. || rando(&colvar.stoch_sample_seed) < colvar.stoch_sample)
	  hills_add(mtd_data);                                      // add HILL
      } 
      if(ntrh) {
	read_hills(mtd_data,0,hills.first_read);              	// is it time to read_hills?
	hills.first_read = 0;
      }
      if(ntwg)  grid_write_tofile(&bias_grid);                         // write GRID on file
    }
    
    cvw.Vwall=soft_walls_engine(colvar.ss0,cvw.fwall);                // Wall potential
    
    Vext=ext_forces_engine(colvar.ss0,&extpot,fext);              // External potential
    
    cvw.Vwall+=steer_engine(colvar.ss0,cvw.fwall);                // Wall potential
    
    cvw.Vwall+=abmd_engine(colvar.ss0,cvw.fwall);                // Wall potential
    
    if (logical.do_dafed) dafed_engine(colvar.ss0);	       // #### d-AFED
    
    cvw.Vwall+=tamd_engine(colvar.ss0,cvw.fwall);                // TAMD/DAFED temperarily added here
  
  
    if (logical.do_constraint) Vconstr+=constraint_engine(mtd_data->dt);            // shake on the cvs 
  
    
    if(logical.do_steerplan)steerplan_engine();               // steerplan to plan adaptive us and more! 
    
    apply_forces(mtd_data);
    
    if(logical.debug_derivatives) debug_derivatives(mtd_data,ntp);
    
    if(ntp) print_colvar_enercv(mtd_data->time);	        	// dump COLVAR

  }

  if(colvar.pg.nlist!=0)calc_projections( &(colvar.pg)); 
   

  stopwhen_engine(); // kill the run in specific points

  if(firstTime)firstTime = 0;			                // the first PluMed step 

}

//------------------------------------------------------------------------------------------

real PREFIX soft_walls_engine(real* ss0,real* force){
  // Wall potential: V_wall(s) = sigma*((s-s0+offset)/redux)**n
  // WARNING: n must be even, sigma>0, redux>0; offset>0; s0 can be positive or negative
  real uscale,lscale,uexp,lexp;
  real V;
  int i;
  V=0.;
  for(i=0;i<colvar.nconst;i++){
    if(force) force[i]=0.0;
    if(logical.upper[i]) {                                                      // if there is a soft wall on this cv
      uscale = (ss0[i]-cvw.upper[i]+cvw.uoff[i])/cvw.ueps[i];                   // calculates the position on the wall 
      if(uscale>0.) {
        uexp = (real) cvw.uexp[i];
        V+=cvw.sigma[i]*pow(uscale, cvw.uexp[i]);
        if(force) force[i]+=-(cvw.sigma[i]/cvw.ueps[i])*uexp*pow(uscale, cvw.uexp[i]-1);
      }
    }
    if(logical.lower[i]) {                                                      // if there is a soft wall on this cv
      lscale = (ss0[i]-cvw.lower[i]-cvw.loff[i])/cvw.leps[i];                   // calculates the position on the wall
      if(lscale<0.) {
        lexp = (real) cvw.lexp[i];
        V+=cvw.lsigma[i]*pow(lscale, cvw.lexp[i]);
        if(force) force[i]+=-(cvw.lsigma[i]/cvw.leps[i])*lexp*pow(lscale, cvw.lexp[i]-1);
      }
    }
  };
  return V;
};
  
//-------------------------------------------------------------------------------------------

real PREFIX tamd_engine(real* ss0,real* force)
{
  real V,tmp;
  int i_c;
  V=0.0;
  if(logical.tamd) for(i_c=0;i_c<colvar.nconst;i_c++)if(colvar.on[i_c]){
    real tmp;
    real ff0,ff1;
    if(firstTime){
// Dummy variable is initialized at the physical temperature
      real c2=sqrt(mtd_data.boltz*tamd.starttemp/tamd.spring[i_c]);
      tmp=c2*rando_gaussian(&tamd.seed);
      tamd.pos[i_c]=ss0[i_c] + tmp;
      fprintf(mtd_data.fplog,"|- TAMD/DAFED %d CV : STARTVALUE %f\n",i_c+1,tamd.pos[i_c]);
      fflush(mtd_data.fplog);
      ff0=ff1=tamd.spring[i_c]*tmp;
    } else {
      tmp=-ss0[i_c]+tamd.pos[i_c];
      if(colvar.type_s[i_c]==5 || ( colvar.type_s[i_c]==34 && colvar.type[i_c]==2 )) {
                       if(tmp > M_PI)
                         tmp -= 2.*M_PI;
                       if(tmp < -M_PI)
                        tmp += 2.*M_PI;
      } 
      real h0,h1;
      ff0=tamd.spring[i_c]*tmp;
      h0=0.5*tamd.spring[i_c]*tmp*tmp;
      real c1=exp(-mtd_data.dt/tamd.tau);
      real c2=sqrt(mtd_data.boltz*tamd.wtemp/tamd.spring[i_c])*(1.0-c1*c1);
      tmp=c1*tmp + c2*rando_gaussian(&tamd.seed);
      ff1=tamd.spring[i_c]*tmp;
      h1=0.5*tamd.spring[i_c]*tmp*tmp;
      tamd.pos[i_c]=ss0[i_c] + tmp;
      tamd.drift+=(h1-h0)*(1/tamd.simtemp-1/tamd.wtemp)/mtd_data.boltz;
    }
    real ff=0.5*(ff1+ff0);
    V+=0.5*tamd.spring[i_c]*tmp*tmp;
    if(force) force[i_c] += ff;
  }
  return V;
};

real PREFIX steer_engine(real* ss0,real* force)
{ 
  real V, tmp;
  real dpos;	 // ### increment in steering reference position.
  real ff=0.;
  int i_c;
#ifdef STANDALONE
  FILE *file;
  char *str, stringa[800], filename[100];
#endif

  V=0.0;
  dpos=0.0;     // ### initializing properly is important for dpos

  for(i_c=0;i_c<colvar.nconst;i_c++) if(logical.steer[i_c]){

    if(firstTime){
      if(cvsteer.impose_start[i_c] == 0) cvsteer.pos[i_c] = cvsteer.start[i_c] = ss0[i_c];
      else cvsteer.pos[i_c] = cvsteer.start[i_c];

      fprintf(mtd_data.fplog,"|- STEERING %d CV : STARTVALUE %f\n",i_c+1,cvsteer.start[i_c]);
      fflush(mtd_data.fplog);
 
      if(cvsteer.max[i_c] < cvsteer.start[i_c]) cvsteer.sign[i_c] = -1; 
      else cvsteer.sign[i_c] = +1;

      // check if you're there since the beginning
      if(cvsteer.pos[i_c]==cvsteer.max[i_c]) {
        fprintf(mtd_data.fplog,"|- STEERING %d CV ARRIVED TO TARGET POINT %f in %d STEPS\n",i_c+1,cvsteer.max[i_c],colvar.it);
        fflush(mtd_data.fplog);
      } 

      // #### Initialize work
      cvsteer.work[i_c]= 0.0;
      cvsteer.old_force[i_c]= 0.0;

#ifdef STANDALONE
      sprintf(filename, "STEER.%d.rst",i_c);
      file = fopen(filename,"w");
      fprintf(file,"%lf  %lf",cvsteer.start[i_c],cvsteer.pos[i_c]);
      fclose(file);
#endif
    } else { // increase when you're not at the first step 
#ifdef STANDALONE 
      // open the file 
      sprintf(filename, "STEER.%d.rst",i_c); 
      file = fopen(filename,"r");
      if(file==NULL){
        char buf[1024];
        sprintf(buf,"Cannot read %s  : EXITING\n",filename);
        plumed_error(buf);
      } else {
        str = fgets(stringa, 800, file); 
        sscanf(str, "%lf %lf",&cvsteer.start[i_c],&cvsteer.pos[i_c]);  
        fprintf(mtd_data.fplog,"|- STEERING %d CV RESTARTED FROM POINT: %f  STARTED: %f\n",i_c+1,cvsteer.pos[i_c],cvsteer.start[i_c]);  
        fclose(file);
        if(cvsteer.max[i_c] < cvsteer.start[i_c]){
          cvsteer.sign[i_c] = -1; 
        } else {
          cvsteer.sign[i_c] = +1;
        } 
      }
      if((logical.not_same_step) && fabs(cvsteer.pos[i_c]-cvsteer.start[i_c])<fabs(cvsteer.max[i_c]-cvsteer.start[i_c])) {
        // #### Introduce variable dpos used for work
        dpos =    cvsteer.sign[i_c] * cvsteer.delta[i_c] / 1000.0;
        cvsteer.pos[i_c] += dpos;
      }
      sprintf(filename, "STEER.%d.rst",i_c); 
      file = fopen(filename,"w");
      fprintf(file,"%lf  %lf",cvsteer.start[i_c],cvsteer.pos[i_c]);
      fclose(file); 
#else
      if((logical.not_same_step) && fabs(cvsteer.pos[i_c]-cvsteer.start[i_c])<fabs(cvsteer.max[i_c]-cvsteer.start[i_c])) {
        // #### Introduce variable dpos used for work
	dpos =  cvsteer.sign[i_c] * cvsteer.delta[i_c] / 1000.0;
    	cvsteer.pos[i_c] += dpos ;
      }
#endif
    } 
    if(fabs(cvsteer.pos[i_c]-cvsteer.start[i_c])>fabs(cvsteer.max[i_c]-cvsteer.start[i_c])) {
      cvsteer.pos[i_c] = cvsteer.max[i_c];
      fprintf(mtd_data.fplog,"|- STEERING %d CV ARRIVED TO TARGET POINT %f in %d STEPS\n",i_c+1,cvsteer.max[i_c],colvar.it);  
      fflush(mtd_data.fplog); 
    }
    /* HERE PUT THE PERIODICITY YOU NEED!!!!!!!! */ 
    tmp = ss0[i_c]-cvsteer.pos[i_c];
    if(colvar.type_s[i_c]==5 || ( colvar.type_s[i_c]==34 && colvar.type[i_c]==2 )) {
      if(tmp > M_PI) tmp -= 2.*M_PI;
      if(tmp < -M_PI) tmp += 2.*M_PI;
    }
    if(cvsteer.annealing[i_c]>0.) { 
      ff = -(cvsteer.annealing[i_c]/mtd_data.temp_t)*cvsteer.spring[i_c]*tmp-(cvsteer.annealing[i_c]/mtd_data.temp_t)*cvsteer.slope[i_c];
      V += 0.5*(cvsteer.annealing[i_c]/mtd_data.temp_t)*cvsteer.spring[i_c]*tmp*tmp+cvsteer.slope[i_c]*(cvsteer.annealing[i_c]/mtd_data.temp_t)*tmp;
    } else {
      ff = -cvsteer.spring[i_c]*tmp - cvsteer.slope[i_c];
      V += 0.5*cvsteer.spring[i_c]*tmp*tmp + cvsteer.slope[i_c]*tmp;
    }
    if(force) force[i_c] += ff;

    // #### Work integration with trapeze scheme
    cvsteer.work[i_c] += 0.5 * dpos * (ff + cvsteer.old_force[i_c]);
    // Alternatively, one could use :
    // work += 0.5 * k * dpos * [  ss0 + ss0_old  - 2*pos + dpos  ]
    cvsteer.old_force[i_c] = ff;
  }
  return V;
}

//---------------------------------------------------------------------------------------------

real PREFIX abmd_engine(real* ss0, real* force)
{
  real V, ff;
  int i_c;

  V=0.0;
  for(i_c=0;i_c<colvar.nconst;i_c++)if(logical.abmd[i_c]){

    abmd.now[i_c] = (ss0[i_c]-abmd.exp[i_c])*(ss0[i_c]-abmd.exp[i_c]);

    if(firstTime&&!logical.restart_abmd) abmd.min[i_c] = abmd.now[i_c];

    if(abmd.now[i_c]<abmd.min[i_c]) abmd.min[i_c] = abmd.now[i_c];
    else {
      ff = -2.*abmd.spring[i_c]*(abmd.now[i_c]-abmd.min[i_c])*(ss0[i_c]-abmd.exp[i_c]);
      if(force) force[i_c] += ff;
      V += 0.5*abmd.spring[i_c]*(abmd.now[i_c]-abmd.min[i_c])*(abmd.now[i_c]-abmd.min[i_c]);
    }

  }
  return V;
}

//---------------------------------------------------------------------------------------------

real PREFIX ext_forces_engine(real* ss0, struct grid_s *grid, real* force)
{
  real Vex = 0.0;

  if(logical.do_external) Vex=grid_getstuff(grid,ss0,force);

  return Vex;
}

//------------------------------------------------------------------------------

void PREFIX init_print_colvar_enercv()
{
/*
  In this routine we print the headers for the COLVAR file.
*/
  int i_c, nactive;
  FILE* cv_file;
  cv_file = fopen((mtd_data.ionode?mtd_data.colfilen:"/dev/null"), "a");
  fprintf(cv_file, "%s", "#! FIELDS");
  fprintf(cv_file, "%s", " time");
  for(i_c=0;i_c<colvar.nconst;i_c++) fprintf(cv_file, " cv%i", i_c+1);
#ifdef CVS
  if( colvar.nbespoke > 0 ){ fprintf(cv_file, " bespoke_error"); }
#endif
  if (logical.do_hills) fprintf(cv_file, " vbias");
  if (logical.do_walls) fprintf(cv_file, " vwall");
  if (logical.do_constraint) fprintf(cv_file, " vconstr");
  if (logical.do_steerplan)  fprintf(cv_file, " vstp");
  if (logical.do_external) fprintf(cv_file, " vext");
#ifdef RECONMETAD
#ifndef DRIVER
//  if( reconOn==1 && reconinpt.monitor!=1 ){ fprintf(cv_file, " vrecon"); }
//  if( reconOn==1 && reconinpt.monitor==1 ){ fprintf(cv_file, " current_basin"); }
  if( reconOn==1 ){ fprintf(cv_file, " vrecon"); }  
#endif
#endif
  for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.steer[i_c]){fprintf(cv_file," XX XX RST%d WORK%d",i_c+1,i_c+1);} }
  for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.abmd[i_c]){fprintf(cv_file," XX XX ABMD%d ",i_c+1);} }
  for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.cnstr[i_c]){fprintf(cv_file," XX XX CONSTRPOS%d CONSTRENE%d ",i_c+1,i_c+1);} }
  if (logical.do_steerplan){ fprintf(cv_file,"%s",steerplan.log)   ;}; 
  if (colvar.pg.nlist!=0)  fprintf(cv_file, "%s",colvar.pg.log ); 
  for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.dafed[i_c]){
    fprintf(cv_file," dAFED%dS dAFED%dT dAFED%dE dAFED%dW ",i_c+1,i_c+1,i_c+1,i_c+1);} }
  for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.tamd && colvar.on[i_c]){fprintf(cv_file," TAMD%d ",i_c+1);} }
  if(logical.tamd) {fprintf(cv_file," TAMD_drift ");};
// end of headers
  fprintf(cv_file,"\n");
// list active CVs (useful e.g. for bias-exchange post-processing)
  if(strlen(colvar.hills_label)>0){
    nactive=0;
    for(i_c=0;i_c<colvar.nconst;i_c++){ if(colvar.on[i_c]) nactive++; }
    fprintf(cv_file, "#! ACTIVE %d",nactive);
    if(nactive>0) {
      for(i_c=0;i_c<colvar.nconst;i_c++){ if(colvar.on[i_c]) fprintf(cv_file, " %d",i_c+1); }
    }
    fprintf(cv_file, " %s",colvar.hills_label); 
    fprintf(cv_file,"\n");
  }
// close COLVAR file
  fclose(cv_file);
}

void PREFIX print_colvar_enercv(real time_s)
{
    int i, i_c; // real consRecon=0.;
    static FILE *cv_file=NULL;

// This will allow us to build a conserved quantity for reconnaissance metadynamics
// #ifdef RECONMETAD
//    if( reconOn==1 ) consRecon=getCons_recon( myreconObj );
// #endif

//#ifdef RECONMETAD 
//#ifndef DRIVER
//    int current_basin=0; 
//    if( reconOn==1 && reconinpt.monitor==1 ){
//      double welikedupes_s[reconinpt.nconst];
//      for(i=0;i<reconinpt.nconst;i++) welikedupes_s[i]=colvar.ss0[reconinpt.cvlist[i]];
//      current_basin=dometa_monitor( reconinpt.nconst, welikedupes_s, myreconObj ); 
//    } 
//#endif
//#endif

    if(!cv_file) cv_file = fopen((mtd_data.ionode?mtd_data.colfilen:"/dev/null"), "a");
/*
ATTENTION: all the quantities written here should be consistent with the header written in init_print_colvar_enercv()
In this way the new colvar format will work properly
To keep the new fmt equivalent to the old one (at least for a transition time)
we leave here also the "comments" such as "RST 3". These numbers are labeled with XX
in the header, and at some point will be removed.
*/
    fprintf(cv_file, "%10.4f", time_s);
    for(i_c=0;i_c<colvar.nconst;i_c++) fprintf(cv_file, "   %14.9f", colvar.ss0[i_c]);
#ifdef CVS
    real bespoke_err=0.; //double cv_in[colvar.bespoke_ncv];
    if( colvar.nbespoke > 0 ){ 
//       for(i=0;i<colvar.bespoke_ncv;i++){ cv_in[i]=colvar.ss0[ colvar.bespoke_cvlist[i] ]; }
       bespoke_err=calculate_bespoke_error( mybespokeObj );
       fprintf(cv_file, "   %14.9f", bespoke_err); 
    }
#endif
    if (logical.do_hills) fprintf(cv_file, " %14.9f ",hills.Vhills/mtd_data.eunit);
    if (logical.do_walls) fprintf(cv_file, " %14.9f ",cvw.Vwall/mtd_data.eunit);
    if (logical.do_constraint) fprintf(cv_file, "  %14.9f ",Vconstr/mtd_data.eunit);
    if (logical.do_steerplan) fprintf(cv_file, "  %14.9f ",Vsteerplan/mtd_data.eunit);
    if (logical.do_external) fprintf(cv_file, "  %14.9f ",Vext/mtd_data.eunit);
#ifdef RECONMETAD
#ifndef DRIVER
    if( reconOn==1 ){ fprintf(cv_file, "   %14.9f", Vrecon/mtd_data.eunit ); }
//    if( reconOn==1 && reconinpt.monitor!=1 ){ fprintf(cv_file, "   %14.9f", Vrecon/mtd_data.eunit ); }
//    if( reconOn==1 && reconinpt.monitor==1 ){ fprintf(cv_file, "   %14d", current_basin ); }
#endif
#endif
    for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.steer[i_c]){fprintf(cv_file," RST %d %14.9f %20.9lf",i_c+1,cvsteer.pos[i_c], cvsteer.work[i_c] );} }
    for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.abmd[i_c]){fprintf(cv_file," ABMD %d %14.9f ",i_c+1, abmd.min[i_c] );} }
    for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.cnstr[i_c]){fprintf(cv_file," CONSTR %d %14.9f %20.9lf ",i_c+1,cvcnstr.pos[i_c],cvcnstr.energy[i_c]);} }
    if(logical.do_steerplan)fprintf(cv_file,"%s",steerplan.log)  ;
    if(colvar.pg.nlist!=0)  fprintf(cv_file, "%s",colvar.pg.log ); 
    for(i_c=0;i_c<colvar.nconst;i_c++){ if(logical.dafed[i_c]){ print_dafed(&dafed[i_c], cv_file, i_c); } }
    for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.tamd && colvar.on[i_c]){fprintf(cv_file," %14.9f ",tamd.pos[i_c]);} }
    if(logical.tamd) {fprintf(cv_file," %14.9f ",tamd.drift);};
    fprintf(cv_file, "\n");
    fflush(cv_file);
#ifdef STANDALONE 
    fclose(cv_file);
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------

void PREFIX commit_analysis()
{
  int i, a, b, ix;
  FILE *commit_file;

  if(firstTime){
    commit_file = fopen((mtd_data.ionode?"COMMIT":"/dev/null"), "a");
    for(i=0;i<colvar.nconst;i++) fprintf(commit_file, "   %14.7f", colvar.ss0[i]);
    fclose(commit_file);
  }

  a = 0;
  b = 0;
  for(i=0;i<commit.ncv;i++){
    ix=commit.index[i];
    if(commit.Amin[ix]<colvar.ss0[ix] && colvar.ss0[ix]<commit.Amax[ix]) a++;
    if(commit.Bmin[ix]<colvar.ss0[ix] && colvar.ss0[ix]<commit.Bmax[ix]) b++;
  }

  if(a==commit.ncv || b==commit.ncv) {
    fprintf(mtd_data.fplog, "|- SYSTEM HAS REACHED AN ENDING REGION\n");
    logical.commit = 0;
    commit_file = fopen((mtd_data.ionode?"COMMIT":"/dev/null"), "a");
    if(a==commit.ncv) fprintf(commit_file, " A \n");
    if(b==commit.ncv) fprintf(commit_file, " B \n");
    fclose(commit_file);
#ifdef PLUMED_GROMACS45
    gmx_set_stop_condition(1); 
#else
    EXIT();
#endif
  }
}

//----------------------------------------------------------------------------------------------------------------------

void PREFIX zero_forces(struct mtd_data_s *mtd_data ) {
  // set to zero all forces
#ifndef PLUMED_GROMACS
  int i;
  for(i=0;i<mtd_data->natoms;i++){
    mtd_data->force[i][0] = 0.0; 
    mtd_data->force[i][1] = 0.0;
    mtd_data->force[i][2] = 0.0;
  }  
#endif

}

void PREFIX apply_forces(struct mtd_data_s *mtd_data )
{
  real ddr, uscale, lscale, dsdt, nor, fact;
  int ix, i, iat, wall,i_c;
  #ifdef NAMD
  Vector f;   
  #endif 

#ifdef RECONMETAD
#ifndef DRIVER
   // Setup an array to contain the derivatives from reconnaissance metadynamics
   real recon_der[colvar.nconst]; for(i=0;i<colvar.nconst;i++) recon_der[i]=0.0;
   if( reconOn==1 ){                             // && reconinpt.monitor!=1 ){
      double welikedupes_s[reconinpt.nconst], welikedupes_d[reconinpt.nconst];
      // Transfer colvars from plumed arrays to dummy array for reconnaissance
      for(i=0;i<reconinpt.nconst;i++) welikedupes_s[i]=colvar.ss0[reconinpt.cvlist[i]];
      // Do the reconnaissance calculation
      real ene; ene=dometa_recon(colvar.it, reconinpt.nconst, welikedupes_s, welikedupes_d, myreconObj);
      // Transfer energy and derivatives to plumed arrays
      Vrecon+=ene; for(i=0;i<reconinpt.nconst;i++) recon_der[reconinpt.cvlist[i]]=welikedupes_d[i]; 
   } 
#endif
#endif

  for(i_c=0;i_c<colvar.nconst;i_c++){

    ddr = colvar.ff_hills[i_c] + cvw.fwall[i_c] + fext[i_c]	// hills, soft wall and external contribution
          +   (real) dafed[i_c].f;               // #### d-AFED contribution

#ifdef RECONMETAD
#ifndef DRIVER
//    if( reconOn==1 && reconinpt.monitor!=1 ) ddr+= recon_der[i_c];      // Reconnaissance metadynamics forces
    if( reconOn==1 ) ddr+=recon_der[i_c];      // Reconnaissance metadynamics forces
#endif
#endif

    if(logical.cnstr[i_c]) ddr-=cvcnstr.lambdadt2[i_c];         // constraint

#ifdef PLUMED_GROMACS
    if(colvar.type_s[i_c]==35) {colvar.d_0[i_c]=ddr; continue;}  // Don't apply force here, do it where you have dd information
#endif

    for(i=0;i<colvar.natoms[i_c];i++) {
      iat = colvar.cvatoms[i_c][i];
#ifdef NAMD 
       f.x=ddr*colvar.myder[i_c][i][0];
       f.y=ddr*colvar.myder[i_c][i][1];
       f.z=ddr*colvar.myder[i_c][i][2];
       addForce(iat,f);
#else
             mtd_data->force[iat][0] += ddr*colvar.myder[i_c][i][0];     // PluMeD forces
             mtd_data->force[iat][1] += ddr*colvar.myder[i_c][i][1];
             mtd_data->force[iat][2] += ddr*colvar.myder[i_c][i][2];
#endif    
    }
  }
}
//
// this engine SHAKEs  the needed d.o.f. (one by one, not interwined)  
// based on Frenkel's book
//
real PREFIX constraint_engine(real tstep0){
  real tstep;
  int go,i,i_c,iat,niter,j;
  real ***posc,***newposc,***velc,***oldposc,***startder;
  FILE *fp;
  char *str,string[100]; 
  real **hack_pos1,totene;
  real lambdadt2,tmp,totlambdadt2;
  char buf[200]; 
  // check whether constraining else skip


  go=0;
  for(i_c=0;i_c<colvar.nconst;i_c++){
     if(logical.cnstr[i_c]){go=1;}
  }
  if(go==0){return 0.;}

  totene=0.;

   //sprintf(buf,"|- ENTERING CONSTRAINING NST %d DT %f \n",mtd_data.istep,tstep );
   //plumed_warn(buf);
# ifdef AMBER
   tstep=tstep0; 
#else
   tstep=tstep0; 
#endif


  // first step: do the mallocs ( always in STANDALONE ) 

#ifndef STANDALONE
  if(mtd_data.istep==0){
#endif

    cvcnstr.newposc  =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.oldposc  =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.velc     =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.posc     =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.startder =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.go       =(int *)malloc(colvar.nconst*sizeof(int));
    cvcnstr.oldcv    =(real *)malloc(colvar.nconst*sizeof(real));

// amber uses its own force vector

// standalone needs alloc
# if defined (STANDALONE)
    cvcnstr.oldforce =(real *)malloc(mtd_data.natoms*3*sizeof(real));
#endif

    for(i_c=0;i_c<colvar.nconst;i_c++){
       if(logical.cnstr[i_c]==1){
         cvcnstr.posc[i_c]       =float_2d_array_alloc(mtd_data.natoms,3); 
         cvcnstr.newposc[i_c]    =float_2d_array_alloc(mtd_data.natoms,3); 
         cvcnstr.oldposc[i_c]    =float_2d_array_alloc(mtd_data.natoms,3); 
         cvcnstr.velc[i_c]       =float_2d_array_alloc(mtd_data.natoms,3); 
         cvcnstr.startder[i_c]   =float_2d_array_alloc(mtd_data.natoms,3); 
       }   
    } 

#ifndef STANDALONE
  }
#endif
  // alias those names 
  posc   =cvcnstr.posc;
  newposc=cvcnstr.newposc;
  oldposc=cvcnstr.oldposc;
  velc=cvcnstr.velc;
  startder=cvcnstr.startder;
  // assign positions
  for(i_c=0;i_c<colvar.nconst;i_c++){
        //fprintf(mtd_data.fplog,"PP %d \n",colvar.natoms[i_c]);  
       if(logical.cnstr[i_c]==1){
          cvcnstr.go[i_c]=1;
          cvcnstr.oldcv[i_c]=colvar.ss0[i_c];
          for(i=0;i<colvar.natoms[i_c];i++){
               iat = colvar.cvatoms[i_c][i];
               posc[i_c][iat][0]=mtd_data.pos[iat][0]; 
               posc[i_c][iat][1]=mtd_data.pos[iat][1]; 
               posc[i_c][iat][2]=mtd_data.pos[iat][2]; 
               startder[i_c][i][0]=colvar.myder[i_c][i][0]; 
               startder[i_c][i][1]=colvar.myder[i_c][i][1]; 
               startder[i_c][i][2]=colvar.myder[i_c][i][2]; 
          }
       }
  }
  if(mtd_data.istep==0){
  // first step: put velocity to zero, put old position to the actual one 
    for(i_c=0;i_c<colvar.nconst;i_c++){
       if(logical.cnstr[i_c]==1){
              for(i=0;i<colvar.natoms[i_c];i++){
                   iat = colvar.cvatoms[i_c][i];
                   oldposc[i_c][iat][0]=posc[i_c][iat][0]; 
                   oldposc[i_c][iat][1]=posc[i_c][iat][1]; 
                   oldposc[i_c][iat][2]=posc[i_c][iat][2]; 
              }
       }
     }
  }
   
#ifdef STANDALONE
  // if standalone check if velocity is available from yout software
  fp = fopen("f.xyz","r");
//  sprintf(buf,"|- TAKING FORCE FILE \n");
//  plumed_warn(buf);
 
  if( fp ) {
    // exists
    i=0;
    while(1){
       str=fgets(string,100,fp);
       if(str==NULL)break;
       if(feof(fp))break;
   //    plumed_warn(str);
       str = strtok( string, " \t" );
       cvcnstr.oldforce[i*3+0]=atof(str);
       str=strtok(NULL," \t");
       cvcnstr.oldforce[i*3+1]=atof(str);
       str=strtok(NULL," \t");
       cvcnstr.oldforce[i*3+2]=atof(str);
       str=strtok(NULL," \t");
       i++; 
    } 
    if(i!=mtd_data.natoms){ plumed_error("|- FORCE FILE DOES NOT MATCH THE NUMBER OF ATOMS ");}
    fclose(fp);
  } else {
    plumed_error("|- NEED FOR FORCE FILE  \n");
  } 
//  sprintf(buf," TAKING VEL FILE \n");
//  plumed_warn(buf);
 
  fp = fopen("v.xyz","r");
  if( fp ) {
    // exists
    for(i_c=0;i_c<colvar.nconst;i_c++){
       if(logical.cnstr[i_c]==1){
          i=0;
          while(1){
             str=fgets(string,100,fp);
             if(str==NULL)break;
             if(feof(fp))break;
          //   plumed_warn(str);
             str = strtok( string, " \t" );
             cvcnstr.velc[i_c][i][0]=atof(str);
             str=strtok(NULL," \t");
             cvcnstr.velc[i_c][i][0]=atof(str);
             str=strtok(NULL," \t");
             cvcnstr.velc[i_c][i][0]=atof(str);
             str=strtok(NULL," \t");
             i++;
          } 
          if(i!=mtd_data.natoms){ plumed_error("|- VELOCITY FILE DOES NOT MATCH THE NUMBER OF ATOMS \n");}
          rewind(fp);
       } 
    }
    fclose(fp);
  } else {
    // doesnt exist ! old positions?
//    sprintf(buf," TAKING OLDPOS FILE \n");
//    plumed_warn(buf);
 
    fp = fopen("oldx.xyz","r");
    if( fp ) {
      for(i_c=0;i_c<colvar.nconst;i_c++){
          if(logical.cnstr[i_c]==1){
             i=0;
             while(1){
                str=fgets(string,100,fp);
                if(str==NULL)break;
                if(feof(fp))break;
             //   plumed_warn(str);
                str = strtok( string, " \t" );
                oldposc[i_c][i][0]=atof(str);
                str=strtok(NULL," \t");
                oldposc[i_c][i][0]=atof(str);
                str=strtok(NULL," \t");
                oldposc[i_c][i][0]=atof(str);
                str=strtok(NULL," \t");
                i++;
             } 
             if(i!=mtd_data.natoms){ plumed_error("|- VELOCITY FILE DOES NOT MATCH THE NUMBER OF ATOMS \n");}
          }
       }
       fclose(fp);
       for(i_c=0;i_c<colvar.nconst;i_c++){
           if(logical.cnstr[i_c]==1){
                  for(i=0;i<colvar.natoms[i_c];i++){
                       iat = colvar.cvatoms[i_c][i];
                       velc[i_c][iat][0]=(posc[i_c][iat][0]-oldposc[i_c][iat][0])/tstep; 
                       velc[i_c][iat][1]=(posc[i_c][iat][1]-oldposc[i_c][iat][1])/tstep; 
                       velc[i_c][iat][2]=(posc[i_c][iat][2]-oldposc[i_c][iat][2])/tstep; 

                       //fprintf(mtd_data.fplog,"|- %d V %f %f %f\n",i,velc[i_c][iat][0],velc[i_c][iat][1],velc[i_c][iat][2]);
                       //fprintf(mtd_data.fplog,"|- %d P %f %f %f\n",i,posc[i_c][iat][0],posc[i_c][iat][1],posc[i_c][iat][2]);
                       // save the old constraint derivative: pay attention on the derivative ordering
                  }
           }
       } 
    } else {
       plumed_error("|- CONSTRAINT NOT POSSIBLE WITHOUT EITHER OLD POSITION OR ACTUAL VELOCITIES \n");
    }
  }
#else
  // calculate velocity through finite difference (take this as default)
  for(i_c=0;i_c<colvar.nconst;i_c++){
      if(logical.cnstr[i_c]==1){
             for(i=0;i<colvar.natoms[i_c];i++){
                  iat = colvar.cvatoms[i_c][i];
// take all the time the findiff so to be sure
//#ifdef AMBER
	// take velocity from the right vector
        //          velc[i_c][iat][0]=cvcnstr.oldvel[iat*3+0]; 
        //          velc[i_c][iat][1]=cvcnstr.oldvel[iat*3+1]; 
        //          velc[i_c][iat][2]=cvcnstr.oldvel[iat*3+2]; 
//#else
	// take velocity from finite differences 
                  velc[i_c][iat][0]=(posc[i_c][iat][0]-oldposc[i_c][iat][0])/tstep; 
                  velc[i_c][iat][1]=(posc[i_c][iat][1]-oldposc[i_c][iat][1])/tstep; 
                  velc[i_c][iat][2]=(posc[i_c][iat][2]-oldposc[i_c][iat][2])/tstep; 
//#endif
             }
      }
  }  

#endif
 
  for(i_c=0;i_c<colvar.nconst;i_c++){
     if(logical.cnstr[i_c]==1){

           if(cvcnstr.verbose[i_c]==1){sprintf(buf," DOING CONSTRAINTS ON CV %d INIT VAL: %f ",i_c+1,cvcnstr.oldcv[i_c]);
           plumed_warn(buf);}
           niter=0;
           totlambdadt2=0.;
           while(cvcnstr.go[i_c]){// check on threshold value
              tmp=cvcnstr.spring[i_c]*(cvcnstr.oldcv[i_c]-cvcnstr.pos[i_c]); // fake constraint spring gradient 
              for(i=0;i<colvar.natoms[i_c];i++){
                  iat = colvar.cvatoms[i_c][i];
                  // make a velocity verlet free evolution (assume the units are congruent)
                  //  r_n= r_o +v dt - 0.5 * dt**2 * F
                  for (j=0;j<3;j++){
#if defined (AMBER)  || defined (STANDALONE)
// time evolution for known codes 
                     newposc[i_c][iat][j]=posc[i_c][iat][j]+  velc[i_c][iat][j]*tstep+ 
                           0.5 * tstep *tstep * 
                     (cvcnstr.oldforce[iat*3+j])  / mtd_data.mass[iat]       
                     - totlambdadt2*tmp*startder[i_c][i][j]/mtd_data.mass[iat];

// I do not know this code
#else
                     sprintf(buf," CONSTRAINTS NOT IMPLEMENTED IN THIS PROGRAM ");
                     plumed_error(buf);
#endif            
                  }
                  
              }
              //calculate the constraint value
              
              // do an hack : unlock the address from mtd_data.pos and redirict to an identical vector newposc[i_c]    

              hack_pos1=mtd_data.pos; // remember this address
              mtd_data.pos=newposc[i_c]; //calculate now with this coordinates 
              switch(colvar.type_s[i_c]){
                  // geometric CVs
                  case 1: dist_restraint(i_c, &mtd_data); break;			// DISTANCE
                  case 2: mindist_restraint(i_c, &mtd_data); break;               	// MINDIST
                  case 3: coord_restraint(i_c, &mtd_data); break;	              	// COORD
                  case 4: angle_restraint(i_c, &mtd_data); break;	                // ANGLE
                  case 5: torsion_restraint(i_c, &mtd_data); break;                  // TORSION
                  case 6: alfabeta_restraint(i_c, &mtd_data); break;                	// ALPHA-BETA
                  // interaction CVs
                  case 7: hbonds_restraint(i_c, &mtd_data); break;                   // HBONDS
                  case 8: dipole_restraint(i_c, &mtd_data); break;       		// DIPOLE
                  // conformations CVs
                  case 11: radgyr_restraint(i_c, &mtd_data); break;   	       	// RGYR
                  case 16: dihcor_restraint(i_c, &mtd_data); break;                 	// DIHEDRAL-COR
                  // water CVs
                  case 20: waterbridge_restraint(i_c, &mtd_data); break;            	// WATERBRIDGE
                  case 21: density_restraint(i_c, &mtd_data); break;            	// WATERBRIDGE
                  case 22: densityswitch_restraint(i_c, &mtd_data); break;           // DENSITYSWITCH
                  // trajectory CVs
                  case 30: spath_restraint(i_c, &mtd_data); break;                   // S_MAPPATH
                  case 31: zpath_restraint(i_c, &mtd_data); break;                   // Z_MAPPATH
                  case 32: position_restraint(i_c, &mtd_data); break;                // ATOM POSITION
                  case 33: elstpot_restraint(i_c, &mtd_data); break;                 // ELSTPOT POSITION
                  case 34: puckering_restraint(i_c, &mtd_data); break;               // PUCKERING
                  case 35: energy_restraint(i_c, &mtd_data); break;                  // ENERGY
                  case 36: helix_restraint(i_c, &mtd_data); break;                   // HELIX
                  case 37: alpharmsd_restraint(i_c, &mtd_data); break;               // ALPHARMSD
                  case 38: antibetarmsd_restraint(i_c, &mtd_data); break;            // ANTIBETARMSD
                  case 39: parabetarmsd_restraint(i_c, &mtd_data); break;            // PARABETARMSD
                  //case 40: camshift_restraint(i_c, &mtd_data); break;              // CAMSHIFT ENERGY
#ifdef CVS
                  case 46: bespoke_restraint(i_c, &mtd_data); break;                 // BESPOKE COLLECTIVE COORDINATES
#endif
                  case 47: rdf_restraint(i_c, 1, &mtd_data); break;                  // DISCRETIZED RDF
                  case 49: histogram_restraint(i_c, &mtd_data); break;               // HISTOGRAM CV
                  case 52: adf_restraint(i_c, 1, &mtd_data); break;                  // DISCRITIZED ANGLE DISTRIBUTION FUNCTION
                  case 55: sprint_restraint(i_c, &mtd_data); break;                   // SPRINT TOPOLOGICAL CV
              }

              mtd_data.pos=hack_pos1;  
              //return     
              //calculate lambda*deltat*deltat 
              tmp=0.;
              for(i=0;i<colvar.natoms[i_c];i++){
                 iat = colvar.cvatoms[i_c][i];
                 tmp+=startder[i_c][i][0]*colvar.myder[i_c][i][0]/mtd_data.mass[iat];
                 tmp+=startder[i_c][i][1]*colvar.myder[i_c][i][1]/mtd_data.mass[iat];
                 tmp+=startder[i_c][i][2]*colvar.myder[i_c][i][2]/mtd_data.mass[iat];
              }   
              // calculate new lambda
              lambdadt2=(colvar.ss0[i_c]-cvcnstr.pos[i_c])/(cvcnstr.oldcv[i_c]-cvcnstr.pos[i_c]);
              lambdadt2/=2*cvcnstr.spring[i_c];
              lambdadt2/=tmp;
              // die?
              tmp=sqrt((colvar.ss0[i_c]-cvcnstr.pos[i_c])*(colvar.ss0[i_c]-cvcnstr.pos[i_c]));  
              if(cvcnstr.verbose[i_c]==1) fprintf(mtd_data.fplog,"|- totlambdadt2 %12.6f lambdadt2 %12.6f tmp %12.6f POS %12.6f CV %12.6f NITER %4d STEP %8.3f\n",totlambdadt2,lambdadt2,tmp,cvcnstr.pos[i_c],colvar.ss0[i_c],niter,tstep);

              //
              // KILL THE ITERATION:
              //
              if(tmp<cvcnstr.delta[i_c])cvcnstr.go[i_c]=0;
              if(niter>cvcnstr.maxiter[i_c])cvcnstr.go[i_c]=0;  

              // apply it  (put here because it allows for a simpler debug when compared with my guineapig program) 
              //if ( cvcnstr.go[i_c]!=0 )totlambdadt2+=lambdadt2;
              totlambdadt2+=lambdadt2;

              niter++;
           }    
           cvcnstr.lambdadt2[i_c]=totlambdadt2*2.*cvcnstr.spring[i_c]*(cvcnstr.oldcv[i_c]-cvcnstr.pos[i_c])/(tstep*tstep); 
           if(cvcnstr.verbose[i_c]==1){sprintf(buf," DONE CONSTRAINT ON CV %d LAMBDA %f ITER %d ",i_c+1,cvcnstr.lambdadt2[i_c],niter);
           cvcnstr.energy[i_c]=totlambdadt2*pow(cvcnstr.oldcv[i_c]-cvcnstr.pos[i_c],2)/(tstep*tstep); 
           totene+=cvcnstr.energy[i_c];
           plumed_warn(buf);}
     }
  }
  // put the actual position as the old position
  for(i_c=0;i_c<colvar.nconst;i_c++){
     if(logical.cnstr[i_c]==1){
        colvar.ss0[i_c]=cvcnstr.oldcv[i_c];
        for(i=0;i<colvar.natoms[i_c];i++){
             iat = colvar.cvatoms[i_c][i];
             oldposc[i_c][iat][0]=posc[i_c][iat][0]; 
             oldposc[i_c][iat][1]=posc[i_c][iat][1]; 
             oldposc[i_c][iat][2]=posc[i_c][iat][2]; 
             colvar.myder[i_c][i][0]=startder[i_c][i][0];
             colvar.myder[i_c][i][1]=startder[i_c][i][1];
             colvar.myder[i_c][i][2]=startder[i_c][i][2];
        }
     }
  }
  //sprintf(buf,"|- EXITING CONSTRAINING \n");
  //plumed_warn(buf);
 
  return totene;
};
 
void PREFIX steerplan_engine()
{ 
   int i,j,k,l,m,thisstage,nextstage,mycv; 
   real time,x1,x2,k1,k2,deltat,tmp,force;
   time=mtd_data.time;
   //fprintf(mtd_data.fplog,"|- STEERPLAN_ON TIME %f FIRSTTIME %d\n",time,firstTime);

   //char string2[200],string[200];
   //for(i=0;i<steerplan.totstages;i++){
   //    fprintf(mtd_data.fplog,"|- TIME: %12.6f  ",steerplan.actions[i].t);  
   //    for(j=0;j<steerplan.ncvs;j++){
   //      if(steerplan.actions[i].activecv[j].wildcardpos){
   //       sprintf(string,"*");
   //      }else {
   //       sprintf(string,"%12.6f",steerplan.actions[i].activecv[j].pos);
   //      }
   //      if(steerplan.actions[i].activecv[j].wildcardk){
   //       sprintf(string2,"*");
   //      }else {
   //       sprintf(string2,"%12.6f",steerplan.actions[i].activecv[j].k);
   //      }
   //      if(steerplan.actions[i].activecv[j].k<0){
   //           fprintf(mtd_data.fplog," CV %3d SKIPPING THIS STAGE... ",(steerplan.actions[i].activecv[j].ncv+1));
   //      }else{
   //           fprintf(mtd_data.fplog," CV %3d KAPPA %12s POS %12s ",steerplan.actions[i].activecv[j].ncv+1,string2,string);
   //      }
   //    }
   //    fprintf(mtd_data.fplog,"\n");  
   //}
 
   // firsttime: define the initial stage and the time for the nextone 
   // which stage I am? New one? calculate new params
#ifdef STANDALONE 
        FILE *file;
        char filename[100],sstr[100] ;
        char *str, stringa[800];
        char buf[1024];
        // open the file 
        if(!firstTime){
          sprintf(filename, "STEERPLAN.rst"); 
          file = fopen(filename,"r");
          if(file==NULL){
            sprintf(buf,"Cannot read %s  : EXITING\n",filename);
            plumed_error(buf);
          }else{
            str = fgets(sstr, 100, file); 
            sscanf(str,"%lf",&steerplan.nextstage_time);
            for(i=0;i<steerplan.ncvs;i++){
                 str = fgets(sstr, 100, file); 
                 sscanf(sstr,"%lf %lf %lf %lf %d",&steerplan.actualcv[i].k,&steerplan.actualcv[i].kv,&steerplan.actualcv[i].x0,&steerplan.actualcv[i].v,&steerplan.actualcv[i].type);
            }
            fclose(file);
          }
        }
#endif
   if(time>=steerplan.nextstage_time||firstTime){
       for(i=0;i<steerplan.totstages;i++){
           if(steerplan.actions[i].t>=time){steerplan.current_stage=i;
              if((i+1)!=steerplan.totstages){steerplan.nextstage_time=steerplan.actions[i+1].t;}
              else{steerplan.nextstage_time=1.e9;} 
           break;} 
       }
       fprintf(mtd_data.fplog,"|- CURRENT STAGE %d  \n",steerplan.current_stage);
       // setup the run
       thisstage=steerplan.current_stage;
       for(i=0;i<steerplan.ncvs;i++){
          mycv=steerplan.actions[thisstage].activecv[i].ncv;
          // is this a skipped stage for this cv?
          if (steerplan.actions[thisstage].activecv[i].k<0){// bypass and keep the values 
                 fprintf(mtd_data.fplog,"|- STEERPLAN CV  %d :NOTHING TO UPDATE TO AT THIS STAGE. I'LL GO ON AS BEFORE\n",steerplan.actions[thisstage].activecv[i].ncv+1);
               // if firstime you don't specify any action this is taken as k=0 pos=* type=1 v=0.  action
               // case1: coming  from nowhere: default is put k=0 (no effective potential)
              if(firstTime){
                       steerplan.actualcv[i].k=0.;
                       steerplan.actualcv[i].v=0.;
                       steerplan.actualcv[i].kv=0.;
                       steerplan.actualcv[i].x0=colvar.ss0[mycv];
                       steerplan.actualcv[i].type=1;
                       // keeep nowhere
              }
          }else{ // scan to the next stages to find the following pinpoints if any
                 fprintf(mtd_data.fplog,"|- STEERPLAN CV  %d SCANNING FOR THE NEXT ACTION\n",steerplan.actions[thisstage].activecv[i].ncv+1);
                 // case2: going to  nowhere: default is keep the old restraint (static: v=0.)
                 // this happens when the list is over or when all the k are <0  
                 if(thisstage==steerplan.totstages-1){
             //    fprintf(mtd_data.fplog,"1 THISSTAGE %d TOTSTAGES %d\n",thisstage,steerplan.totstages);
                       if(!steerplan.actions[thisstage].activecv[i].wildcardk){steerplan.actualcv[i].k=steerplan.actions[thisstage].activecv[i].k;}
                       steerplan.actualcv[i].v=0.;
                       steerplan.actualcv[i].kv=0.;
                       if(steerplan.actions[thisstage].activecv[i].wildcardpos){
                                steerplan.actualcv[i].x0=colvar.ss0[mycv]; 
                       }
                       else{ steerplan.actualcv[i].x0=steerplan.actions[thisstage].activecv[i].pos;}
                       steerplan.actualcv[i].type=steerplan.actions[thisstage].activecv[i].type;
                 }else{ 
                       for(j=thisstage+1;j<steerplan.totstages;j++){
                    //       fprintf(mtd_data.fplog,"XXCV %d K %f WC %d\n",mycv+1,steerplan.actions[j].activecv[i].k,steerplan.actions[j].activecv[i].wildcardk);
                           if((steerplan.actions[j].activecv[i].k>=0.) || (steerplan.actions[j].activecv[i].wildcardk==1)){
                           break;} // got it!!!!
                       } 
                       if(j==steerplan.totstages){//no more actions to take! 
                           if(!steerplan.actions[thisstage].activecv[i].wildcardk){steerplan.actualcv[i].k=steerplan.actions[thisstage].activecv[i].k;}
                           steerplan.actualcv[i].v=0.;
                           steerplan.actualcv[i].kv=0.;
                           if(steerplan.actions[thisstage].activecv[i].wildcardpos){
                               steerplan.actualcv[i].x0=colvar.ss0[mycv];
                           }
                           else{ steerplan.actualcv[i].x0=steerplan.actions[thisstage].activecv[i].pos;}
                           steerplan.actualcv[i].type=steerplan.actions[thisstage].activecv[i].type;
                    //       fprintf(mtd_data.fplog,"2 THISSTAGE %d TOTSTAGES %d\n",thisstage,steerplan.totstages);
                       }else{ 
                    //      fprintf(mtd_data.fplog,"3 THISSTAGE %d TOTSTAGES %d\n",thisstage,steerplan.totstages);
                          deltat=(steerplan.actions[j].t-steerplan.actions[thisstage].t)/mtd_data.dt;              
                          // copy the actual center (wildcard?) : start from the actual real cv 
                          if(steerplan.actions[thisstage].activecv[i].wildcardpos){steerplan.actualcv[i].x0=colvar.ss0[mycv];}
                          else{ steerplan.actualcv[i].x0=steerplan.actions[thisstage].activecv[i].pos;}
                          // copy the next center (wildcard?) : start from the previous center!!!
                          if(steerplan.actions[j].activecv[i].wildcardpos){
                              steerplan.actualcv[i].v=0.;   
                              steerplan.actualcv[i].nowhere=1;
                          }
                          else { 
                             x1=steerplan.actualcv[i].x0; 
                             x2=steerplan.actions[j].activecv[i].pos  ;
                             steerplan.actualcv[i].v=(x2-x1)/deltat   ;
                             steerplan.actualcv[i].nowhere=0;
                          };
 
                           // copy the actual k (wildcard? dont do anything ) 
                          if(!steerplan.actions[thisstage].activecv[i].wildcardk){steerplan.actualcv[i].k=steerplan.actions[thisstage].activecv[i].k;}
                          // copy the next k (wildcard?) : keep this one 
                          if(steerplan.actions[j].activecv[i].wildcardk){k2=steerplan.actions[thisstage].activecv[i].k;}
                          else { k2=steerplan.actions[j].activecv[i].k  ;};
                          k1=steerplan.actualcv[i].k;
                          steerplan.actualcv[i].kv=(k2-k1)/deltat   ;
                          steerplan.actualcv[i].type=steerplan.actions[thisstage].activecv[i].type;
                     //     fprintf(mtd_data.fplog,"XXX CV %d X1 %f X2 %f K1 %f K2 %f \n",steerplan.actions[thisstage].activecv[i].ncv+1,x1,x2,k1,k2); 
                       } 
                 } 
                  
          }
       }
   }else{// evolove normally
      for(i=0;i<steerplan.ncvs;i++){
            steerplan.actualcv[i].k=steerplan.actualcv[i].k+steerplan.actualcv[i].kv;
            steerplan.actualcv[i].x0=steerplan.actualcv[i].x0+steerplan.actualcv[i].v;
            if(steerplan.actualcv[i].nowhere){
                     mycv=steerplan.actions[steerplan.current_stage].activecv[i].ncv;
                     steerplan.actualcv[i].x0=colvar.ss0[mycv]; 
            }

      }      
   }
#ifdef STANDALONE
       sprintf(filename,"STEERPLAN.rst");
       file = fopen(filename,"w");
       fprintf(file,"%lf\n",steerplan.nextstage_time);
       for(i=0;i<steerplan.ncvs;i++){
            fprintf(file,"%lf %lf %lf %lf %d\n",steerplan.actualcv[i].k,steerplan.actualcv[i].kv,steerplan.actualcv[i].x0,steerplan.actualcv[i].v,steerplan.actualcv[i].type);
       }
       fclose(file);
#endif

   // create a log string
   strcpy(steerplan.log," STP ");
   char cvlog[100]; 
   for(i=0;i<steerplan.ncvs;i++){
          mycv=steerplan.actions[steerplan.current_stage].activecv[i].ncv;
          //fprintf(mtd_data.fplog,"|- STEERPLAN CV  %d : \n",steerplan.actions[steerplan.current_stage].activecv[i].ncv+1);
          //fprintf(mtd_data.fplog,"|- X  %12.6f XV %12.6f  K %12.6f KV %12.6f  \n",steerplan.actualcv[i].x0,steerplan.actualcv[i].v,steerplan.actualcv[i].k,steerplan.actualcv[i].kv);
          sprintf(cvlog," CV %2d X %12.6f K %12.6f T %1d ",mycv+1,steerplan.actualcv[i].x0,steerplan.actualcv[i].k,steerplan.actualcv[i].type);
          strcat(steerplan.log,cvlog);
   }  
   //fprintf(mtd_data.fplog,"%s\n",steerplan.log);

    // now the plan is fine: do what you need!!

   for(i=0;i<steerplan.ncvs;i++){
      mycv=steerplan.actions[steerplan.current_stage].activecv[i].ncv; 
      tmp=colvar.ss0[mycv]-steerplan.actualcv[i].x0;
        
      if(steerplan.actualcv[i].type==2){
         // positive: wall only if positive
         if(tmp<0.)tmp=0.;
      } else if (steerplan.actualcv[i].type==3){
         // negative: wall only if negative 
         if(tmp>0.)tmp=0.;
      } 
      if(colvar.type_s[mycv]==5 || ( colvar.type_s[mycv]==34 && colvar.type[mycv]==2 )){
                       if(tmp > M_PI)
                         tmp -= 2.*M_PI;
                       if(tmp < -M_PI)
                        tmp += 2.*M_PI;
      } 
      //force = -steerplan.actualcv[i].k*tmp;
      steerplan.actualcv[i].force = -steerplan.actualcv[i].k*tmp;
      steerplan.actualcv[i].energy = steerplan.actualcv[i].k*tmp*tmp*0.5;
      cvw.fwall[mycv] +=  steerplan.actualcv[i].force ;
      Vsteerplan  +=  steerplan.actualcv[i].k*tmp*tmp*0.5; 
   } 

//   fprintf(mtd_data.fplog,"|- STEERPLAN_ON END\n");
}

//---------------------------------------------------------------------------------------------
void PREFIX stopwhen_engine (){
    int i;
    int stop;
    stop=0;
    for(i=0;i<colvar.nconst;i++){
      if(stopwhen.actmin[i]){
              if(colvar.ss0[i]<stopwhen.min[i]){
                  fprintf(mtd_data.fplog, "|-PLUMED STOPWHEN:\n" );
                  fprintf(mtd_data.fplog, "|-NOW KILLING BEGAUSE CV  %d IS LESS THAN %f :ACTUAL VALUE IS %f \n",i+1,stopwhen.min[i],colvar.ss0[i] );
                  fflush(mtd_data.fplog);
                  stop=1;
              }
      }
      if(stopwhen.actmax[i]){
              if(colvar.ss0[i]>stopwhen.max[i]){
                  fprintf(mtd_data.fplog, "|-PLUMED STOPWHEN:\n" );
                  fprintf(mtd_data.fplog, "|-NOW KILLING BEGAUSE CV  %d IS MORE THAN %f : ACTUAL VALUE IS %f \n",i+1,stopwhen.max[i],colvar.ss0[i] );
                  fflush(mtd_data.fplog);
                  stop=1;
              }
      }
    } 
    if(stop){
#ifdef PLUMED_GROMACS45
                  gmx_set_stop_condition(1); 
#else
                  EXIT();
#endif
    }
}
//---------------------------------------------------------------------------------------------

