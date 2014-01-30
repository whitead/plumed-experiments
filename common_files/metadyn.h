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

#ifndef EXTERNALS
#if defined NAMD
   #define MYEXT 
#elif defined LAMMPS_PLUMED
   #define MYEXT 
#elif defined PLUMED_CPMD
   #define MYEXT
#else
   #define MYEXT extern
#endif   
#else
   #define MYEXT
#endif

#ifdef RECONMETAD 
  #include "recon_cbind.h"
#endif

// common header files
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#ifdef HAVE_MATHEVAL
#include <matheval.h> 
#endif

// NAMD header files and definitions
#if defined (NAMD)
#define PREFIX GlobalMasterMetaDynamics::
#include "ComputeHomePatches.h"
#include "GlobalMaster.h"
#include "GlobalMasterEasy.h"
#include "NamdTypes.h"
#include "SimParameters.h"
#include "Molecule.h"
#include "Node.h"
#include "Vector.h" 
#include "signal.h"

// GROMACS header files and definitions
#elif defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
#define PLUMED_GROMACS
#define PREFIX
#include "config.h"
#include "typedefs.h"
#include "smalloc.h"
#include "pbc.h"
#include "vec.h"
#include "physics.h"
#include "network.h"
#include "nrnb.h"
#include "bondf.h"
#include "random.h"
#include "repl_ex.h"

#ifdef PLUMED_GROMACS45
#include "sighandler.h"
#include "xtcio.h"
#ifdef GMX_LIB_MPI
#define MPI
#endif
#else
#ifdef GMX_MPI
#define MPI
#endif
#endif

//ACEMD header files and definitions
#elif defined (ACEMD)
#define PREFIX
#include <aceplug.h>

//LAMMPS header files and definitions
#elif defined (LAMMPS_PLUMED)
#include "domain.h"
using namespace LAMMPS_NS;
#define PREFIX Plumed::

//CPMD
#elif defined (PLUMED_CPMD)
#define PREFIX Plumed::

// Fortran codes (no header needed)
#else
#define PREFIX
#endif

// Eventually, we include the mpi header
#ifdef MPI
#include "mpi.h"
#endif

// This is needed because of the inconsistency in the definition of gmx_repl_ex_t in versions 3 and 4
#if defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
    typedef gmx_repl_ex_t plumed_repl_ex;
#endif

// here we provide alternatives of gromacs macros for the other codes
#if ! defined (PLUMED_GROMACS)
typedef double real;
typedef double rvec[3];
#define REAL_EPS 1.11022302E-16
// double precision accuracy 
#if defined (PLUMED_AIXXLF)
#define snew(ptr,nelem) (ptr)= (nelem==0 ? NULL : (__typeof__(ptr)) calloc(nelem,sizeof(*(ptr))))
#define srenew(ptr,nelem) (ptr)= (__typeof__(ptr)) realloc(ptr,(nelem)*sizeof(*(ptr)))
#else
#define snew(ptr,nelem) (ptr)= (nelem==0 ? NULL : (typeof(ptr)) calloc(nelem,sizeof(*(ptr))))
#define srenew(ptr,nelem) (ptr)= (typeof(ptr)) realloc(ptr,(nelem)*sizeof(*(ptr)))
#endif
#else
// gromacs "real" precision accuracy 
#define REAL_EPS GMX_REAL_EPS
#endif
// use gromacs and lammps proprietary sfree
#if !defined   (PLUMED_GROMACS) && !defined (LAMMPS_PLUMED)
#define sfree(ptr) if(ptr != NULL)free(ptr)   
#endif


// common global structures and definitions
// JFD>
// Changed by recommendation of Michael McGovern.
#define nconst_max 25					// fixed maximum number of COLVARS
// <JFD
#define DP2CUTOFF 6.25					// sigma^2 considered for gaussian
#define GTAB 1000000					// mesh for exponential tablature
#ifndef M_PI
#define M_PI 3.14159265
#endif
#ifndef M_2PI
#define M_2PI 6.2831853
#endif
// JFD>
// Added to implement the McGovern-de Pablo boundary consistent hills
#define M_sqrt2oPI 0.79788456     // sqrt(2/pi)
#define M_sqrt2 1.41421356        // sqrt(2)
// <JFD
// path dimensions
#define MAXATOMS_PATH 900
#define NMAX_PATH 10 
#define MAXFRAMES_PATH 52 
#define MAXATOMS_RMSD 900
#define MAXCHARS_PATH 40
// cmap 
#define MAXDIM_CMAP 5000
#define MAXNUM_GROUP  10
#define MAXATOM_GROUP 30
// stackdimension for hills
#define STACKDIM  10000 
// verlet list
#define MAXNN 2000

// Structure containing the parsed input file
typedef struct {
  int     nlines;     // number of lines
  int*    nwords;     // number of words (in each line)
  char*** words;      // words in each line (+ an extra empty word as a terminator(
  int     ngroups;    // number of group definitions found
  char**  groupnames; // names of groups
  int*    natoms;     // number of members (for each group)
  int**   atoms;      // members (for each group)
} t_plumed_input;

// structure for gradient projections
// linked list
struct  coupling_ll {
        int *at1,*at2,nat1,nat2;
        struct coupling_ll  *next_elem;
};
// couple of cv containing the linked list
struct  el_couple {
    int cv1,cv2;
    struct coupling_ll *first_elem;
};
struct el_diagonal {
    struct coupling_ll *first_elem;
    real *accm;
    real *acct;
};
struct proj_grad_s {
  int *list;
  int  nlist; 
  int nvar,ncouples;
  struct el_couple *couple;
  struct  el_diagonal *diagonal;
  real  **matrix;
  real  **invmatrix;
  real  volume;
  real  *averages;
  real  *prec;
  char dumpfile[200];
  char log[1000];
  real w_stride;
};
// pdb struct for pdb parsing
struct pdb {
	int natoms,*index, *resid;
	real *x,*y,*z,*occ,*beta;
	char **name, **resname,**chain;
};

// use this for generic rmsd container and general path variables
// all dynamical allocatable for ease of use
struct rmsd_container_t {
	// this keeps the label and the original file
	struct pdb *mypdb;
	int ndisplace;
	int nalign;
	int natoms;
	int dummy;
	int dref;
	int dref_freq;
	char dreffile[200];
	int *index;
	real *x,*y,*z;
	real *align; //[MAXATOMS_PATH];    // for each atom is >0 if used for alignment, else 0
	real *displace; //[MAXATOMS_PATH]; // for each atom is >0 if used for displacement, else 0
	real walign;
	real wdisplace;
	int *lalign;
	int *ldisplace;
	int simple; 
	FILE *fpdreffile;
};

// common NAMD/GROMACS interface
struct mtd_data_s
{
   int  natoms;
   real **pos;
   real **vel;
   real **force;
   real *charge;
   real *mass;
   real temp_t; 
   int  repl;
   int  nrepl;
   int  repl_ex_nst;  
   real rte0;
   real rteio;
   real time;
   real time_offset;
   FILE *fplog;
   char metaFilename[120];
   real dt;
   long long int  istep;					// #### This needs to be long long for long RESPA runs !
   long long int  istep_old;				// #### This needs to be long long for long RESPA runs !
   real eunit;
   real boltz;
   int  ionode;                                        // this is true only in the master node
   real energy;
   int  newcolvarfmt;
   char   colfilen[800];					// COLVAR and ENERCV files
   char   hilfilen[800];					// HILLS file
   char   basefilen[800];				// BASENAME file
   char   dump_filen[800];                               // file for atoms dumping
#if !defined (PLUMED_GROMACS45)
   int    dump_file;                                     // unit
#else
   t_fileio* dump_file;
#endif
   int    dump_stride;                                   // stride
   int    dump_atoms;                                   // list of dumped atoms
   int*   dump_list;                                    // list of dumped atoms
#ifdef MPI
// communicators for parallel PLUMED
   MPI_Comm comm;       // INTRA replica
   MPI_Comm intercomm;  // INTER replica (for bias-exchange and ptmetad)
#endif

// code-specific definitions
#if defined (PLUMED_GROMACS)
   const t_commrec *mcr;  
   real cell[3][3];
   t_pbc  metapbc;
#if defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
   int  ePBC;
#endif

#elif defined (OPEP)
   char log[120];
   int imcon;

#elif defined (DL_POLY) || defined (GAT_LJONES)
   int imcon;
   real cell[9];

#elif defined (ACEMD) || defined (AMBER) || defined (DRIVER) || defined (STANDALONE) || defined (PLUMED_QESPRESSO)
   int  imcon;
   real cell[3];

#elif defined (LAMMPS_PLUMED)
   Domain *mydomain;
#endif
#ifdef STANDALONE 
   real ampli;
#endif

  
};

struct logical_s
{
  int    do_hills;					// hills on/off
  int    widthadapt;
  int    do_inversion;                                  // do inversion condition 
  int    invert[nconst_max][2];                         // inversion condition
  int    restart_hills;					// restart meta 
  int    append;                                        // append on COLVAR 
  int    restart_abmd;                                 	// restart abmd
  int    upper[nconst_max];				// upper walls
  int    lower[nconst_max];				// lower walls
  int    interval[nconst_max];                          // fahimeh 
  int    do_walls;						// #### any walls ?
  int    steer[nconst_max];				// steering cv on/off
  int    dafed[nconst_max];				// d-AFED cv on/off ####
  int    do_dafed;					// d_AFED main switch ###
  int    do_constraint;					// constraint main switch ###
  int    tamd;                                          // tamd/dafed on/off
  int    abmd[nconst_max];				// abmd cv on/off
  int    cnstr[nconst_max];				// constraint cv on/off
  int    always[nconst_max];				// if one in (hills, walls, steering, abmd, constraint) is on then is 1 else is 0
  int    remd;						// replica exchange parallel tempering
  int    rpxm;						// replica exchange metadynamics
  int    hrex;                                          // Hamiltonian replica-exchange
  int    commit;					// committors analysis
  int    print;						// do COLVAR 
  int    norescale;
  int    welltemp;
  // JFD>
  int    dicksonian_tempering;      // Use a well-tempered metadynamics such that Dickson's approximation is exact.
  int    transition_tempering;      // Use transition-tempering
  int    ttdebug;
  int    mcgdp_hills;          // Use the McGovern-de Pablo boundary consistent hills.
  // <JFD
  //ADW>
  int    target_distribution; //do target distribution metadynamisc
  // <ADW
  int    lreflect[nconst_max];
  int    ureflect[nconst_max];
  int    debug;
  int    not_same_step;
  int    meta_inp;
  int    parallel_hills;                                // hills sum is parallelized
  int    debug_derivatives;
  int    enable_untested_features;
  int    do_grid;
  int    read_grid;
  int    write_grid;
  int    donot_spline;
  int    debug_grid;
  int    do_walkers;
  int    puckering;
  int    path;
  int    energy;
  int    read_old_bf;
  int    do_external;
  int    do_alphabetarmsd;
  int    do_sprint;
  int    do_steerplan;
  int    do_pca;
  int    nlist[nconst_max];
};

struct adapt
{
  int    block;
  int    stride;
  int    on;
  real   inf;
  real   sup;
  real   widthmultiplier;                               // Hills width adapt
};

struct colvar_s
{
  int nbespoke;                                      // # of Bespoke CVs
#ifdef CVS
  // int    bespoke_eps;                                   // The exponent used to calculate the bespoke cvs
  // int    bespoke_nneigh;                                // The number of neighbours used to calculate bespoke cvs
  int    bespoke_ncv;                                   // # number of CVs in each bespoke CV (these should all be the same)
  int    *bespoke_cvlist;                               // The CVs used to calculate each bespoke collective coordinate                         
#endif
  int    nconst;					// # CVs
  int    nt_print;					// step for printing colvar and enercv
  int    it;						// PluMeD step counter
  real   delta_r [nconst_max];				// Hills width along CVs
  real   inv_limit[nconst_max][2];                      // Limits used in inversion
  real   inv_ref[nconst_max];                           // Reflection interval used in the inversion
  real   inv_inv[nconst_max];                           // Inversion interval
  real   inv_maxww[nconst_max];                         // Inversion gaussian height upper limit (colvar.inv_maxww*hills.wwr)
  real   **delta_s;					// Hills width along CVs in time
  struct adapt adapt[nconst_max];                             // adaptive width structure
  real   ff_hills[nconst_max];				// force due to hills
  int    on      [nconst_max];				// hills on/off on a CV
  int    type_s  [nconst_max];				// colvar type (DIST, ...)
  // Parameters for histogram colvars
  int    histo_ncv[nconst_max];                         // The number of cvs used to make each histogram
  int    *histo_cvlist[nconst_max];                     // The list of cvs used to make the histogram
  real   histo_low[nconst_max];                         // The lower bound of the histogram bead  
  real   histo_high[nconst_max];                        // The upper bound of the histogram bead
  real   histo_width[nconst_max];                       // The widths of the Gaussians used to make histograms
  // Parameters for rdf colvars
  int    rdfNorm [nconst_max];                          // Are we normalizing according to the volume of the bin?
  int    rdflab  [nconst_max];                          // Which rdf is this bead from
  real   rdfBeadLower[nconst_max];                      // Lower bound for bin 
  real   rdfBeadUpper[nconst_max];                      // Upper bound for bin
  real   rdfBeadWidth[nconst_max];                      // Width of bin
  // End of rdf colvar params
  int    doTrig  [nconst_max];                          // Used by torsion so that we can use sines and cosines of torsions GAT 
  int    nn      [nconst_max];				// used by for numerator exponent
  int    mm      [nconst_max];				// used by for denominator exponent
  real   r_0     [nconst_max];				// used by for binding distance
  real   d_0     [nconst_max];				// used by for binding distance
  real   beta    [nconst_max];				// used by mindist
  int    intpar  [nconst_max][21];                      // array of integers (general use)
  rvec   realpar [nconst_max][21];                      // array of reals (general use)
  rvec   vecpar  [nconst_max][21];                       // array of vecors (general use)
  real   *map0   [nconst_max]; 				// inter e intra contact starting maps
  real   *map1   [nconst_max];                          // inter e intra contact starting maps
  int    groups  [nconst_max];				// energy groups id, other id
  int    type    [nconst_max];				// type id (beta sheet, alpha elicas, none, ecc..)
  rvec   *myder  [nconst_max];				// derivatives
  real   ss0     [nconst_max];                          // CVs value 
  real   Mss0    [nconst_max];				// Hills width adapt
  real   M2ss0   [nconst_max];				// Hills width adapt
  real   wtemp;	         				// well tempered temperature
  real   simtemp;                       // simulation temperature
  real   wfactor;                       // welltemp factor = wtemp/simtemp
  //ADW>
  int    b_treat_independent  [nconst_max];          //
  //<ADW
  // JFD>
  real   tttemp;                        // transition tempered temperature
  real   ttfactor;                      // transitiontemp factor = tttemp/simtemp
  real   ttthreshold;                   // transitiontemp threshold
  int    n_transition_wells;
  int    transition_wells_converted;
  real   transition_wells[100][nconst_max];                    // Transition wells in transition tempered metadynamics
  int    transition_well_indices[100];
  char   ttdebug_file[800];
  // <JFD
  int    list[nconst_max][4];                           // structure definition for list 
  int    natoms  [nconst_max];
  int    *cvatoms[nconst_max];
  int    logic[nconst_max];
  int    cell_pbc[nconst_max];                          // switch for applying pbc (where it applies)
  int    ptmetad_neighbours;
  real    ptmetad_sigma;
  int    align_atoms;
  int    *align_list;
  struct rmsd_container_t rmsd_container[nconst_max];  // general rmsd container
  // projections 
   struct proj_grad_s  pg;
  // optimized weights
  real   *ow_weight  [3*nconst_max];
  real   hrex_energy;
  // PCA CV
  int  pca_align_atoms;					// number of atoms to align with (*)
  int  *pca_align_list;					// list of atoms to align with (*)
  rvec *pca_align_coord;				// coordinates of atoms to align with (*)
  real  *pca_align_mass;				// masses of atoms to align with (*)
  real  pca_align_totmass;				// total mass of atoms to align with (*)
							// (*) it is not an array [nconst_max] because only 1 reference
							// frame is allowed even if there are more PCA CVs
  struct pcacomp_s *pcacomp[nconst_max];		// pca vector components structure 
  int pca_diff		[nconst_max];			// DIFF keyword flag
  int pca_align		[nconst_max];			// NOALIGN keyword flag
  int pca_upstride	[nconst_max];			// UPSTRIDE keyword flag
  real old_d[3][3];					// saves the rotation matrix for efficiency
  real old_dd_dr1[3][3][3][MAXATOMS_RMSD];		// and its derivatives!

  char  hills_label[50];                                // a label associated to the hills file, shown in headers
};

struct hills_s
{
  real     wwr;						// Hill height
  real     rate;					// Hill deposition rate
  real     max_height;    				// Maximum height of added hills (0 means NO maximum)
  int      max_stride;    				// Maximum stride between hills (0 means NO maximum)
  real     *ww;						// Hill height history
  long int n_hills;					// Hills added
  int      nt_hills;					// Period in step to add hills
  int      nr_hills;					// Period in step to read hills
  real     **ss0_t;					// Hills center history
  char     dir[800];					// HILLS place
  long int ntothills;					// max number of hills
  real     exp[GTAB];					// table for exponential
  long int read;
  fpos_t   *line_counter;
  int      first_read;
  int      nwalkers;
  int      idwalker;
  real   Vhills;					// Hills potential
  // JFD>
  // Extra information for depositing McGovern-de Pablo hills
  int mcgdp_reshape_flag[nconst_max];
  real hill_upper_bounds[nconst_max];
  real hill_lower_bounds[nconst_max];
  real n_hills_near_lower_bound;
  real n_hills_near_upper_bound;
  real erf[GTAB + 1][nconst_max];
  // <JFD
};

struct wall
{
  real   upper   [nconst_max];				// upper limit where start the wall
  real   lower   [nconst_max];				// lower limit
  real   lsigma  [nconst_max];				// lower force constant
  real   sigma   [nconst_max];				// upper force constant
  real   fwall   [nconst_max];				// force due to wall
  int    uexp    [nconst_max];				// upper softwall exponent
  int    lexp    [nconst_max];				// lower softwall exponent
  real   ueps    [nconst_max];				// redux factor for upper wall 
  real   leps    [nconst_max];				// redux factor for lower wall 
  real   uoff    [nconst_max];                          // offset for upper wall
  real   loff    [nconst_max];                          // offset for lower wall     
  int    st_inv  [nconst_max];
  real   Vwall;						// Wall potential
};

struct interval_s
{
  real   upper_limit   [nconst_max];                          // fahimeh 
  real   lower_limit   [nconst_max];                          // fahimeh    
};

struct steer {
 real    pos     [nconst_max];				// position of umbrella
 real    delta   [nconst_max];				// increment of umbrella along cv
 real    spring  [nconst_max];				// elastic constant of umbrella
 real    max     [nconst_max];				// limit
 real    start   [nconst_max];				// start position
 real    slope   [nconst_max];                          // additional linear potential slope*(s-pos)
 real    annealing[nconst_max];                         // if the temperature changes the constants are rescaled as constant=factor/temp_now
 int     sign    [nconst_max];
 int     impose_start [nconst_max];                      // logical for imposing a starting point
 real 	  old_force [nconst_max];                       // #### Old value of the force for work integration
 double  work   [nconst_max];				// #### Work of the steering potential
} ;

struct tamd {
 real    pos     [nconst_max];                          // position of the restraint
 real    spring  [nconst_max];                          // strength of the restraint (from sigma)
 real    tau;                                           // relaxation time
 real    simtemp;                                       // T
 real    wfactor;                                       // T'/T
 real    wtemp;                                         // T'/T
 int     seed;
 real    starttemp;                                     // temperature for pos initialization
 real    drift;
};

struct constraint {
 real    pos     [nconst_max];				// position of umbrella
 real    delta   [nconst_max];				// tolerance
 real    spring  [nconst_max];				// elastic constant of umbrella
 int     maxiter [nconst_max]; 
 real    force   [nconst_max];				// force due to constraint 
 real    lambdadt2 [nconst_max];				// force due to constraint 
 real    energy [nconst_max];				// force due to constraint 
 int    verbose [nconst_max];				// verbosity 
 // some other stuff connected with coordinate
 real ***posc,***newposc,***velc,***oldposc,***startder;
 int *go;
 real *oldcv; 
 // amber interface for old force
#if defined (AMBER) || defined (STANDALONE)
 real   *oldforce; 
 real   *oldvel; 
#endif
} ;

// verlet list 
struct  nlist_s{
  real rcut;
  real rskin;
//  int  step;
  int  *nn;
  int  **ni;
  real **base;
};

struct abmd {
 real    exp     [nconst_max];                          // ideal destination
 real    spring  [nconst_max];                          // elastic constant 
 real    min     [nconst_max];                          // best value reached
 real    now     [nconst_max];                          // start position
} ;

struct grid_s {
 real       min      [nconst_max];                      // GRID inferior limit
 real       max      [nconst_max];                      // GRID superior limit
 real       lbox     [nconst_max];                      // GRID bin size
 real       oldelta  [nconst_max];                      // store old HILLS delta 
 real       dx       [nconst_max];                      // GRID spacing 
 real       minilbox [nconst_max];                      // redux GRID bin size 
 real       *pot                 ;                      // array for meta bias
 real       **force              ;                      // array for forces
 real       cutoff               ;                      // store genereal DP2CUTOFF
 real       mem                  ;                      // memory info
 int        bin      [nconst_max];                      // number of bins in total GRID
 int        minibin  [nconst_max];                      // number of bins in redux GRID
 int        period   [nconst_max];                      // periodic ?
 int        index    [nconst_max];                      // to map back the id of active CV
 int        size                 ;                      // size of total GRID 
 int        minisize             ;                      // size of redux GRID
 int        **one2multi          ;                      // from 1d index to multidimensional for redux GRID
 int        **one2multi_full     ;                      // same for full GRID
 int        ncv                  ;                      // number of ACTIVE CVs 
 int        w_stride             ;                      // stride for GRID writing on file
 long int   nhills               ;                      // Total number of HILLS put on GRID
 char       r_file[800]          ;                      // GRID file to read from
 char       w_file[800]          ;                      // GRID file to write to
} ;

struct commit_s {
 int    ncv                 ;                          // number of ACTIVE CVs 
 int    index   [nconst_max];                          // to map back the id of active CV
 real   Amin    [nconst_max];                          // A state for committors analysis
 real   Amax    [nconst_max];                          // ""
 real   Bmin    [nconst_max];                          // B ""
 real   Bmax    [nconst_max];                          // ""
} ;

struct stopwhen_s {                // truncate the run and dump the checkpoint (for ffs) only in gromacs  
 int  actmin[nconst_max];   // is it active on the minimum
 int  actmax[nconst_max];   // is it active on the maximum 
 real  min[nconst_max];      // value min  
 real  max[nconst_max];      // value max
};

// d-AFED ########## -----------------------------------------------------------
struct ggmt_s {
 /* This is a second-order Generalized Gaussian Moments Thermostat */
 /* Static parameters */
 double   Q1;			// thermostat mass1
 double   Q2;			// thermostat mass2
 double   kT0;			// reference temperature times the Boltzmann factor
 double   dt[4];		// sub-time steps for internal RESPA and the Suzuki-Yoshida factorization
 int   	   n_respa_ggmt;	// number of iterations for internal RESPA (default 1)
 /* Dynamic variables */
 double		v1;		// velocity1
 double		v2;		// velocity2
 double 	eta1;		// position1
 double		eta2;		// position2
};

// dafed_s is different from steer :
// steer is an structure of arrays
// dafed will be defined as an array of structures
struct dafed_s {
 /* Static parameters */
 double	 temperature ;		// temperature for meta-variable
 double	 mass   ;		// mass of meta-variable
 double kappa   ;		// coupling constant
 double tauthermo  ;		// time constant for meta-thermostat
 int     n_respa;		// Number of RESPA integration steps between two MD steps
 double	 dt_respa;		// sub-time step for DAFED RESPA
 double periodicity_low ;               // lower limit for periodic boundary conditions on s
 double periodicity_high;               // upper limit for periodic boundary conditions on s
 double periodicity_gap;                // upper limit for periodic boundary conditions on s
 /* Switches */
 int	 do_initialize_s;		// if we set initial conditions
 int	 do_skip_integration;   // if we skip the first DAFED integration step after reading restart
 int     do_periodicity;                // if we use periodic boundary conditions on s
 /* Dynamic variables */
 double	 colvar;		// collective variable in double precision.
 double s      ;		// position of meta-variable
 double vs	;		// velocity of meta-variable
 double	 f ;			// DAFED force
 struct ggmt_s	 ggmt	;	// Generalized Gaussian Moments Thermostat
 double	 dafed_work;		// Work of meta_variable on physical system
 double dWold;			// Old work element for trapeze integration
 double	 thermo_work;		// Thermostat work
 int 	 do_jacobian_force;	// Apply force compensating the Jacobian
 double jacobian_force;		// The value of the Jacobian force
};

struct dafed_control_s {
 int   n_respa;				// Number of RESPA integration steps between two MD steps
 int   write_freq;			// write frequency. zero: do not write.
		 					// -1 write at same time as state.cpt files
 int   do_cpt;				// switch to detect when gromacs writes cpt file
 int   restart; 			// read the restart file?
 char  in_file[800];		// filename for reading
 char  out_file[800];		// filename for writing
};

// ########## -----------------------------------------------------------

struct mathfunction_s{
// if you have it not defined then creates a fake struct 
     void *f, **f_prim;    
     char fline[200]; 
     int  *indcvs;    
     real  *valcvs;    
     char **names;
     int  count;
};


struct rmsd_inpack{
       int natoms;
       real r0[3][MAXATOMS_RMSD];
       real r1[3][MAXATOMS_RMSD];
       real mass[MAXATOMS_RMSD];
       real totmass;
};

struct rmsd_outpack{
       real r0p[3][MAXATOMS_RMSD];//centered reference  frame
       real r1p[3][MAXATOMS_RMSD];//centered running frame  
       real cmr0[3]; //center of mass of reference frame
       real cmr1[3]; //center of mass of running frame
       real err;
       real derr_dr0[3][MAXATOMS_RMSD];
       real derr_dr1[3][MAXATOMS_RMSD];
       real dderr_dr1_dr1[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
       real dderr_dr0_dr0[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
       real dderr_dr1_dr0[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
       real dderr_dr0_dr1[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
       real d[3][3];
       real dd_dr0[3][3][3][MAXATOMS_RMSD];
       real dd_dr1[3][3][3][MAXATOMS_RMSD];
};
//
// a dynamic structure for the workarrays
// 
struct rmsd_workstruct_s{
	int maxsize,natoms;
	real **r0,**r1;
	real *align,*displace;
	real walign,wdisplace;
	int *lalign,*ldisplace;
	int nalign,ndisplace,simple;
	real **r0p ; //[3][MAXATOMS_RMSD];//centered reference  frame
	real **r1p ; //[3][MAXATOMS_RMSD];//centered running frame  
	real **r0p_rotated ; // rotated frames
	real **r1p_rotated ; // 
	real cmr0[3]; //center of mass of reference frame
	real cmr1[3]; //center of mass of running frame
	real err;
	real **derr_dr0; //[3][MAXATOMS_RMSD];
	real **derr_dr1; //[3][MAXATOMS_RMSD];
	real d[3][3];
	real dinv[3][3];
	real ****dd_dr0; //[3][3][3][MAXATOMS_RMSD];
	real ****dd_dr1; //[3][3][3][MAXATOMS_RMSD];
	// optional for 2nd der
	int maxsize_secondder;
	real ****dderr_dr1_dr1; //[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
	real ****dderr_dr0_dr0; //[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
	real ****dderr_dr1_dr0; //[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
	real ****dderr_dr0_dr1; //[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
	// work arrays for rmsd 
	real ****dm_r1_store; //[4][4][3][MAXATOMS_RMSD];
	real ****dm_r0_store;// [4][4][3][MAXATOMS_RMSD];
	real **array_3_n; //[3][MAXATOMS_RMSD];
	//real **derr_dr1_tmp; //[3][MAXATOMS_RMSD];
	//real **derr_dr0_tmp; //[3][MAXATOMS_RMSD];
	real ****dd_dr_temp; //[3][3][3][MAXATOMS_RMSD];	
};

struct rmsd_mini_outpack{
     real r0p[3][MAXATOMS_RMSD];//centered reference  frame
     real r1p[3][MAXATOMS_RMSD];//centered running frame  
     real cmr0[3]; //center of mass of reference frame
     real cmr1[3]; //center of mass of running frame
     real err;
     real derr_dr0[3][MAXATOMS_RMSD];
     real derr_dr1[3][MAXATOMS_RMSD];
     real d[3][3];
     real dd_dr0[3][3][3][MAXATOMS_RMSD];
     real dd_dr1[3][3][3][MAXATOMS_RMSD];
};
// For path in contact map space
struct group_struct{
       int number;
       int numatom[MAXNUM_GROUP];
       int index[MAXNUM_GROUP][MAXATOM_GROUP];
       int index_to_list[MAXNUM_GROUP][MAXATOM_GROUP];
       real rcm[MAXNUM_GROUP][3];
};
struct cmap_pack {
       int index1[MAXDIM_CMAP];
       int index2[MAXDIM_CMAP];       
       int index_from1[MAXDIM_CMAP];
       int index_from2[MAXDIM_CMAP];
       int atoms;
       int list[MAXATOMS_PATH];
       int nn[MAXDIM_CMAP];
       int nd[MAXDIM_CMAP];
       int number;
       int gnumber;
       real r0[MAXDIM_CMAP];
       real weight[MAXDIM_CMAP];
       real cutoff[MAXDIM_CMAP];
       real cmap[MAXFRAMES_PATH][MAXDIM_CMAP];
       int logical_group;
       struct group_struct group;
       };
struct cmap_inpack{
       real r0[MAXATOMS_PATH][3];
       real cmap[MAXDIM_CMAP];
};
struct cmap_outpack{
       real err;
       real derr_dr0[3][MAXATOMS_PATH];
       real derr_dcm[MAXDIM_CMAP];
};
struct coordinates_frameset {
        int natoms;
        int resid[MAXATOMS_PATH],atmnum[MAXATOMS_PATH];
        char resname[MAXATOMS_PATH][MAXCHARS_PATH];
        char label[MAXATOMS_PATH][MAXCHARS_PATH];
        real pos[MAXATOMS_PATH][3];
        char residue[MAXATOMS_PATH];
        int frameset_to_coord[MAXATOMS_PATH]; // for each atoms in the frameset provides the corresponding index on running coord
        int frameset_to_align[MAXATOMS_PATH]; // for each atoms in the frameset provides the corresponding index used alignment(-1 if not involved in rmsd) 
        int align_to_coord[MAXATOMS_PATH];
        int align_to_frameset[MAXATOMS_PATH];
        int ndisplace;
        int nalign;
        real align[MAXATOMS_PATH];    // for each atom is >0 if used for alignment, else 0
        real displace[MAXATOMS_PATH]; // for each atom is >0 if used for displacement, else 0
        real walign;
        real wdisplace;
        int simple; 
};
struct jac_index_s{
	int elem1,elem2,cv1,cv2,natom1,natom2,pos1,pos2;
};
// this is the final unti that contains all the single datas
struct hybrid_elem{
		// which cv
        int cvtype,cvindex;
		//
       	// distance: distance is only one float, but requires the full list of atoms so to copy the derivative 
		//
        real    *ref_dist;    // reference dist for this frame: would work for all scalars 
        int		nder;        // ref number of atoms in the derivative 
        int     *backtable;  // ref intdex of the atoms in the derivative (back transfer of forces)         
	    rvec    *myder ;	//    derivatives: contains the number of atoms involved
		int		ncv ; //one single call can calculate many cv (as rmsd)
		int		contains_pdb;
		struct  pdb *pdb;		// dynamic pdb structure
		// 
		rvec	**cvder;	// one cv per variable. Also needed for Jacobian
		int		*natoms_per_cv;
		int		**list_per_cv;
		char    remarkfield[100];
		// only to be allocated if needed
		real    **diff_der1,**diff_der2,**diff_der3;	
		int		has_diff_alloc,has_jac_alloc;
};
struct hybrid_frameset{
       	// this should contain all the required data (possibly dynamically allocated)
       	// for all the structures needed in the hybrid frameset definition
        // one frame may contain many structures
		int    hbd_nelem;        //numbr of calls to cv calculators
		int    hbd_totalatoms;   // total number of atoms involved
		int    hbd_totcvs;		 // total number of cvs
		int    fixed;			 // flag: is it fixed?
        struct hybrid_elem **hbd_elem; // dynamic structure to allocate the element
        int    *backtable;
        rvec   *myder;
		// the tracking structure for the jacobian (only used by the running frame)
		struct jac_index_s *jac_index;
		int	   jac_index_length,has_jac_alloc;
}; 

struct sz_data {
        int    number;
        real lambda;
        struct coordinates_frameset **frameset;
        struct hybrid_frameset **hbd_frameset;
        struct hybrid_frameset  *hbd_running; //assumin that the number of atoms in each reference  stays constant
        char   names[MAXCHARS_PATH];
        int grad_on,umb_on,mass_on,targeted_on,sqrt_on,norot_on,nocenter_on,indexing_type;
        real gradmax,gradmin,gradk;
        int umblagsteps,umbblocksize,umbstride,umbcount,umbpermanency,countperm;
        real umbtolerance;
        real ****umb_block;// dimensions:  3,natoms,nframes,nblock
        real ***umb_avg;// dimensions:  3,natoms,nframes
        real ***umb_map_block;// dimensions:  tot_con,nframes,nblock
        real **umb_map_avg;// dimensions:  tot_con,nframes
        // if you use debug routine for derivative respect to the path-> I want to keep it safe in the struct
        real ***dpath_dr;
        // just allocated whenever preprocesed with  
        struct cmap_pack my_cmap_pack;
        char   path_type[10];
        int nneigh,*lneigh,neigh,neigh_time;
        // hybrid path
        int *lhybrid,nhybrid,*lcvhybrid; 
        real **mathybrid;
		real ***mathybrid_blocks;
		int do_blocks;
		real **myinverse;
        // bernd ensing path on the fly evolution 
        int  iforce; // switch on force and evolution
        int  ievol; // switch on force and evolution
        int  reset; //reseting wcount (weights) to zero
        real fadefactor;       
        real **disp;
        real *wcount;
        char evol_path[100];
        FILE *fp_evol;
		int debug_metrics;
		int intraframe_dist;
		char intraframe_dist_inputfile[100];
		char intraframe_dist_outputfile[100];
		int intraframe_diff;
	    char intraframe_diff_outputfile[100];

};


struct couplingmatrix_s{
	struct sz_data mysz;
	int dumpfreq;
	int is_on;
	char filename[100];
	FILE *fp;
};

// reference distance matrices for CVs alpharmsd, antibetarmsd, parabetarmsd
struct ref_dist_mat_s{
  real alpha[30][30];
  real alpha_pairs;
  real antibeta[30][30];
  real antibeta_pairs;
  real parabeta[2][30][30];
  real parabeta_pairs[2];
};

// data for sprint CVs
struct sprint_data_s{
  int nat;
  int icv;
  int step;
  real *lambda;
  real **cm;
  real ***grad;
};

// steerplan

struct steeronecv_s{
        real k,pos;
        // take from endpoint 
        int  ncv,wildcardpos,wildcardk;
        // type central=0;positive=1;negative=2; 
        int type;
};

struct steeraction_s{
        struct steeronecv_s  *activecv; 
        real t;
};

struct steeractual_s{
       real v,k,x0,kv,force,energy;
       int type,nowhere;
};

// steerplan.action[stage].cv 
struct steerplan_s {
  int current_stage,totstages;
  struct steeraction_s *actions;
  struct steeractual_s *actualcv;
  real nextstage_time;
  int ncvs; 
  int isactive[nconst_max];
  char log[600]; 
}; 

// pca CV
struct pcacomp_s {
        int a1;                 // atom 1 (read, fixed)
        rvec coeff;             // coefficient (read, fixed)
};

// camshift stuff
struct cam_shift_s {
  char cam_data[nconst_max][128];
  char cam_ff[nconst_max][64];
  int  disu[nconst_max];
  int  num[nconst_max];
  int  mumo[nconst_max];
  double grains[nconst_max];
};


// NAMD CLASS DEFINITION AND PECULIAR METHODS
#ifdef NAMD
class GlobalMasterMetaDynamics : public GlobalMasterEasy 
{

 public:
 GlobalMasterMetaDynamics(); 
 virtual void easy_init(const char *);
 void easy_calc(void);
 void mtd_data_init(  );
 void init_metadyn( );
 void rvec2vec(rvec rv,Vector *v);
#elif LAMMPS_PLUMED
class Plumed
{
 public:
 Plumed(char *metainp, char *metaout , int *atoms, real *mss, real *chg, real *dt, real myboltz); //constructor
 void mtd_data_init(char *metainp, char *metaout, int atoms, real *mss, real *chg, real *dt , real myboltz );
 void init_metadyn( char *metainp, char *metaout, int *atoms, real *mss, real *chg, real *dt , real myboltz);
 void meta_force_calculation(int *allidx, rvec *allpos,rvec *allforce, int allnum , Domain *domain);
 void sfree(void *ptr);
#elif defined(PLUMED_GROMACS)
 void mtd_data_init (int ePBC, real *charge, real *mass,
                    int natoms, real dt, int repl_ex_nst, int repl,
                    int nrepl, real rte0, real rteio, const t_commrec *mcr, FILE *fplog);

 void init_metadyn(int natoms, int ePBC, real *charge, real *mass,
                  real dt, int repl_ex_nst, plumed_repl_ex repl_ex,
                  const t_commrec *mcr, FILE *fplog);
 void ptmetad(real *Epota, real *Epotb, real *Epotba, real *Epotab, int a, int b);
 void ptmetad_vbias(int,real*,real*);
 void ptmetad_helper();
 void ptmetad_exchfluct(int);
 void bias_exchange_traj(int nrepl, int *seed, int *ind);
 void meta_force_calculation(int, int, real (*pos)[3], real (*force)[3], real box[3][3], real energy, real temp);
#if defined (PLUMED_GROMACS4)
 void plumed_setstep(long long int istep, bool bCPT);
#elif defined (PLUMED_GROMACS45)
 void plumed_setstep(long long int istep, gmx_bool bCPT);
#else
 void plumed_setstep(int istep);
#endif
 int  plumed_dd_index(int iat);
 int  plumed_ll_index(int iat);
#if defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
#define oprod cprod
#endif
#elif OPEP
 void init_metadyn_(int *atoms, real *ddt, int *pbc_opep,
                   int *repl, int *nrepl,real *rte0, real *rteio, real *mass,
                   char *lpath, char *logfile, char *metainp, int ll, int mm, int jj);
 void meta_force_calculation_(real *pos, real *force);
 real pbc_mic_(real *rin);
 void mtd_data_init(int pbc, real tstep,int atoms, int repl, int nrepl, real rte0, real rteio, real *mass, char *lpath, char *logfile, char *metainp);
 void ptmetad_vbias(int,real*,real*);
 void share_bias_(int *rep, real *Vbias, real *Vbiasx); 
 void bias_exchange_(int *nrep, int *biaseed, int *ind);
 void bias_exchange_traj(int nrepl, int *seed, int *ind);
 void switch_fluct_(int *rep); 
 void ptmetad_exchfluct(int);
#elif ACEMD 
 void mtd_data_init( real *charge, real *mass,
                    int natoms, real dt, int repl,
		     int nrepl, real rte0, real rteio, char *metainp, real *box);
 void init_metadyn(int natoms, real *charge, real *mass, 
                   real dt, int repl, int nrepl, 
                   real rte0, real rteio, char *metainp, real box[3]);
 void meta_force_calculation(struct aceplug_sim_t* );
#elif defined(DL_POLY) || defined(AMBER) || defined(PLUMED_QESPRESSO) || defined (GAT_LJONES)
 void mtd_data_init(int atoms, real dt ,real *mass, real *charge, int *imcon, real *eunit, char *metainp);
 void meta_force_calculation_(real *cell, int *istep, real *xxx, real *yyy, real *zzz, real *fxx, real *fyy, real *fzz, real *energy);
 void init_metadyn_(int *atoms, real *ddt, real *mass, real *charge, int *imcon, real *eunit, char *metain, int pp);
 void images_(int *i,int *j,int *k,int *natoms,real *cell,real *xxx,real *yyy,real *zzz); 
#elif defined(PLUMED_CPMD)
// these two are the wrapper rotines, needed to call c++ code from fortran
#if defined(PLUMED_CPMD_NOUNDERSCORE)
// some IBM-AIX compilers (not all of them) do not want the underscore...
extern "C" void init_metadyn(int *atoms, int *nsp, int *na, int *nsx, int *nax, real *ddt, int *now, real *mass);
extern "C" void meta_force_calculation(real *pos, real *force, int *nsp, int *na, int *nax, int *nsx);
extern "C" void pbc_cpmd_plumed( real*, real*, real* );
#else
extern "C" void init_metadyn_(int *atoms, int *nsp, int *na, int *nsx, int *nax, real *ddt, int *now, real *mass);
extern "C" void meta_force_calculation_(real *pos, real *force, int *nsp, int *na, int *nax, int *nsx);
extern "C" void pbc_cpmd_plumed_( real*, real*, real* );
#endif
// this is the c++ class containing all plumed (like this there is no name conflict between things in cpmd and things in plumed)
class Plumed{
 public:
 Plumed(); //constructor
 void mtd_data_init( int atoms, int nsp, int *na, int nsx, int nax, real ddt, int now, real *mass, char *metainp);
 // these two are inside the class, not visible from outside, called from the wrapper routines above
 void init_metadyn(int *atoms, int *nsp, int *na, int *nsx, int *nax, real *ddt, int *now, real *mass);
 void meta_force_calculation(real *pos, real *force, int *nsp, int *na, int *nax, int *nsx);
//#elif RECON_DRIVER
// void init_metadyn_(int *atoms, real *ddt, real *mass, real *charge, int *pbc, real *box, char *metainp, int* ncv, double *periods, real *w, real *height, real *sizeParam, int ll);
// void mtd_data_init(int atoms, real *mass, real *charge, char *metainp, int pbc, real *box, real ddt);
// void cv_calculation_(real *box, real *pos, int *ncv, real *cv);
// void ptmetad(real *Epota, real *Epotb, real *Epotba, real *Epotab, int a, int b);
// void ptmetad_vbias(int,real*,real*);
// void ptmetad_sharepot(int nrepl, int repl);
// void bias_exchange_traj(int nrepl, int *seed, int *ind);
// void ptmetad_exchflut(int repl);
// void ptmetad_exchfluct(int);
#elif DRIVER
#define min(a, b) (((a) < (b)) ? (a) : (b))
 void init_metadyn_(int *atoms, real *mass, real *charge, int *pbc, real *box, real *ddt ,char *metainp, int* ncv, int ll);
 void mtd_data_init(int atoms, real *mass, real *charge, char *metainp, int pbc, real *box, real ddt);
 void cv_calculation_(real *box, real *pos, int *ncv, real *cv);
 void ptmetad(real *Epota, real *Epotb, real *Epotba, real *Epotab, int a, int b);
 void ptmetad_vbias(int,real*,real*);
 void ptmetad_sharepot(int nrepl, int repl);
 void bias_exchange_traj(int nrepl, int *seed, int *ind);
 void ptmetad_exchflut(int repl);
 void ptmetad_exchfluct(int);
#elif STANDALONE 
#define min(a, b) (((a) < (b)) ? (a) : (b))
 void init_metadyn_(int *atoms, real *mass, real *charge, int *pbc, real *box, real *tstep, int *nstep, real *myboltz, real *ampli, char *metainp, int ll);
 void mtd_data_init(int atoms, real *mass, real *charge, int pbc, real *box,  real *tstep, int *nstep , real* myboltz, real *ampli , char *metainp );
 void cv_calculation_standalone_(real *box, real *pos, real *force, real *ene);
#endif

 real rando_gaussian(int *ig);
#if ! defined (PLUMED_GROMACS)
// Vector operation (locally reimplemented for codes either than GROMACS)
 real rando(int *ig);
 void oprod(const rvec a,const rvec b,rvec c);
 real iprod(const rvec a,const rvec b);
 real norm(const rvec a);
 real norm2(const rvec a);
 real cos_angle(const rvec a,const rvec b);
 real dih_angle(rvec xi, rvec xj, rvec xk, rvec xl,
               rvec r_ij,rvec r_kj,rvec r_kl,rvec m,rvec n,
               real *cos_phi,real *sign);
 void clear_rvec(rvec a);
#endif

// common declarations
#if defined (PLUMED_GROMACS) && defined (__cplusplus)
extern "C" {
#endif
 void EXIT();
 void read_restraint(struct mtd_data_s *mtd_data);
 void read_defaults(  );
 void plumed_read_input(t_plumed_input* input,FILE* file,FILE* log);
 void plumed_clear_input(t_plumed_input*);
 int  plumed_get_group(const char *word,int **atoms,int n,t_plumed_input*,FILE *log);
 void plumed_error(const char*);
 void plumed_warn(const char*);
 int  plumed_atoi(const char* word);
 double plumed_atof(const char* word);
 int  plumed_get_words(char* line,char*** words);
 void plumed_sum    (struct mtd_data_s *mtd_data,int n,real* buffer);
 void plumed_sumi   (struct mtd_data_s *mtd_data,int n,int* buffer);
 void plumed_intersum(struct mtd_data_s *mtd_data, int n, real *buffer);
 int  plumed_comm_size(struct mtd_data_s *mtd_data);
 int  plumed_comm_rank(struct mtd_data_s *mtd_data);
 int  seek_word(char **word, const char *wanted);
 int  seek_word2(char **word, const char *wanted, int is);
// DAFED routines ########################################################
  void initialize_ggmt(struct ggmt_s *ggmt,double T,double tau, double delta_t);
  real energy_ggmt(struct ggmt_s ggmt);
  void initialize_dafed(struct dafed_s *daf, real dt);
  void dafed_engine(real *this_colvar);
  void print_dafed(struct dafed_s *daf, FILE *cv_file, int i_c);
  void write_dafed_state();
  void read_dafed_state();
  void get_kt(struct dafed_s *daf, double *d);
  void integrate_ggmt(struct dafed_s *daf);

#if defined (PLUMED_GROMACS) && defined (__cplusplus)
}
#endif
// PBS
 void minimal_image(rvec pos1, rvec pos2, real *mod_rij, rvec rij);
// HILLS stuff
 void hills_add(struct mtd_data_s *mtd_data);
 void hills_push(struct mtd_data_s *mtd_data,real ww, real* ss,real* delta);
 real hills_engine(real*,real*);
 real hills_engine_dp(int ih,real* ss0,real* dp);
 void hills_force();
 //Sepearated zero forces from apply forces
 //ADW>
 void zero_forces(struct mtd_data_s *mtd_data );
 //<ADW
 void apply_forces(struct mtd_data_s *mtd_data );
 void inversion_on_boundaries(struct mtd_data_s *mtd_data,int ncv);
 real soft_walls_engine(real*,real*);
 real steer_engine(real*,real*);
 real abmd_engine(real*,real*);
 real tamd_engine(real*,real*);
 real ext_forces_engine(real* ss0, struct grid_s *grid, real* force); 
 real constraint_engine(real tstep);
 void stopwhen_engine(); 
 void steerplan_engine();
 void steer_cv(int);
 void init_print_colvar_enercv();
 void print_colvar_enercv(real time_s);
 void hills_adapt();
 void commit_analysis();
 void read_hills(struct mtd_data_s *mtd_data, int restart, int first_read);
 void hills_reallocate(struct mtd_data_s *mtd_data);
 // JFD>
 // Added for transition tempering
 double transition_bias_ND();
 #define TTMETAD_LATTICE_PBC 0
 #define TTMETAD_LATTICE_RBC 1

 // The lattice struct holds the grid dimensions, number of grid points in each
 // dimension, the boundary conditions in each dimension, and calculation intermediates.

 typedef struct {
     size_t dimensions;
     size_t total_size;
     size_t *dim_sizes;
     size_t *dim_strides;
     size_t *dim_bcs;
 } lattice;

 // Lattice array structs hold an array of data defined on each lattice point and
 // a pointer to the lattice structure of the data.

 typedef struct {
     lattice *index_lat;
     double *vals;
 } lattice_array;

 // The indexed_value_max_heap struct stores indexed values in a heap sorted
 // by value so that the root entry has the greatest value in the heap.

 typedef struct {
     size_t max_size;
     size_t curr_size;
     size_t *indices;
     double *vals;
 } indexed_value_max_heap;

 // Pure lattice routines
 lattice *make_lattice(size_t dims, size_t *d_sizes, size_t *d_bcs);
 void free_lattice(lattice *lat);
 size_t *one_d_index_to_multi_d(lattice *lat, size_t index);
 size_t multi_d_index_to_one_d(lattice *lat, size_t *index);
 void find_neighbors(lattice *lat, size_t index, size_t *n_neighs, size_t *neighs);

 // Lattice-point data array routines
 lattice_array *make_lattice_array(lattice *lat);
 lattice_array *make_const_lattice_array(lattice *lat);
 void free_lattice_array(lattice_array *lat_arr);
 void free_const_lattice_array(lattice_array *lat_arr);
 void set_lattice_array(lattice_array *lat_arr, double *vals);
 void set_const_lattice_array(lattice_array *lat_arr, double *vals);

 // Min heap routines for storing indexed values sorted by value
 indexed_value_max_heap *make_indexed_value_max_heap(size_t max_size);
 void free_indexed_value_max_heap(indexed_value_max_heap *ival_mheap);
 void insert_indexed_value(indexed_value_max_heap *ival_mheap, size_t index, double value);
 double get_max(indexed_value_max_heap *ival_mheap);
 size_t get_max_index(indexed_value_max_heap *ival_mheap);
 void delete_max(indexed_value_max_heap *ival_mheap);
 void swap_heap_entries(indexed_value_max_heap *ival_mheap, size_t index_one, size_t index_two);

 // A routine to find the highest minimum value of a given positive function
 // (passed tabulated in a lattice_array) along any paths between a source
 // and a sink on the lattice.
 double find_maximal_path_minimum(lattice_array *lat_arr, size_t source, size_t sink);
 // <JFD
// CV routines
 //ADW>
 void independent_insert_hack(int i_c, int atom_index);
 void independent_remove_hack(int i_c, int atom_index);
 //ADW<
 void restraint(struct mtd_data_s *mtd_data);
 void test_derivatives(struct mtd_data_s *mtd_data);
 void poly_restraint_testder(int i_c, struct mtd_data_s *mtd_data);
 void func_restraint_testder(int i_c, struct mtd_data_s *mtd_data);
#if defined RECONMETAD
 void test_recon_derivatives(struct mtd_data_s *mtd_data);
#endif
#if defined CVS
 void bespoke_test_derivatives(struct mtd_data_s *mtd_data); 
#endif
 void histogram_testderivatives(int j_c, struct mtd_data_s *mtd_data);
 void debug_derivatives(struct mtd_data_s *mtd_data,int);
// reading...
 int  read_dist          (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_mindist       (char **word,int count,t_plumed_input *input,           FILE *fplog); 
 int  read_coord         (char **word,int count,t_plumed_input *input,           FILE *fplog); 
 int  read_angle         (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_torsion       (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_hbonds        (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_dipole        (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_rgyr          (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_waterbridge   (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_density       (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_densityswitch (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_path          (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_position      (char **word,int count,t_plumed_input *input,           FILE *fplog); 
 int  read_elstpot       (char **word,int count,t_plumed_input *input,           FILE *fplog); 
 int  read_puckering     (char **word,int count,t_plumed_input *input,           FILE *fplog); 
 int  read_energy        (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_alpharmsd     (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_antibetarmsd  (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_parabetarmsd  (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_cmap          (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_pca           (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_sprint        (char **word,int count,t_plumed_input *input,           FILE *fplog);
#ifdef CVS
 int  read_bespoke       (char **word,int count,t_plumed_input *input,           FILE *fplog);
#endif 
 int  read_histogram     (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_rdf           (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_adf           (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_msd           (char **word,int count,t_plumed_input *input,           FILE *fplog);
// these variables read lines directly from the input, thus they need the iline pointer to increment it
 int  read_alfabeta      (char **word,int count,t_plumed_input *input,int *iline,FILE *fplog);
 int  read_dihcor        (char **word,int count,t_plumed_input *input,int *iline,FILE *fplog);
 int  read_helix         (char **word,int count,t_plumed_input *input,int *iline,FILE *fplog);
 int  read_steerplan     (char **word,int count,t_plumed_input *input,int *iline,FILE *fplog);
 int  read_poly          (char **word,int count,t_plumed_input *input,int *iline,FILE *fplog);
 int  read_func          (char **word,int count,t_plumed_input *input,int *iline, int *nw ,FILE *fplog);


 int  load_steerplan     (char *word,struct steerplan_s *steerplan,FILE *fplog);
 void read_couplingmatrix  ( char **word, t_plumed_input *input,FILE *fplog );
	void calc_couplingmatrix  ( int i );

// calculating
#ifdef CVS
 void bespoke_restraint       (int i_c, struct mtd_data_s *mtd_data);
#endif
 void rdf_restraint           (int i_c, int ignore_repeats, struct mtd_data_s *mtd_data);
 void adf_restraint           (int i_c, int ignore_repeats, struct mtd_data_s *mtd_data);
 void histogram_restraint     (int i_c, struct mtd_data_s *mtd_data);
 void dist_restraint          (int i_c, struct mtd_data_s *mtd_data);
 void pt_from_axis_restraint  (int i_c, struct mtd_data_s *mtd_data);
 void pt_from_axis_inplane_restraint (int ptype, int i_c, struct mtd_data_s *mtd_data);
 void proj_on_axis_restraint  (int i_c, struct mtd_data_s *mtd_data);
 void diffdist_restraint      (int i_c, struct mtd_data_s *mtd_data);
 void mindist_restraint       (int i_c, struct mtd_data_s *mtd_data);
 void coord_restraint         (int i_c, struct mtd_data_s *mtd_data);
 void angle_restraint         (int i_c, struct mtd_data_s *mtd_data);
 void torsion_restraint       (int i_c, struct mtd_data_s *mtd_data);
 void alfabeta_restraint      (int i_c, struct mtd_data_s *mtd_data);
 void hbonds_restraint        (int i_c, struct mtd_data_s *mtd_data);
 void dipole_restraint        (int i_c, struct mtd_data_s *mtd_data);
 void radgyr_restraint        (int i_c, struct mtd_data_s *mtd_data);
 void dihcor_restraint        (int i_c, struct mtd_data_s *mtd_data);
 void waterbridge_restraint   (int i_c, struct mtd_data_s *mtd_data);
 void density_restraint       (int i_c, struct mtd_data_s *mtd_data);
 void densityswitch_restraint (int i_c, struct mtd_data_s *mtd_data);
 void spath_restraint         (int i_c, struct mtd_data_s *mtd_data);
 void zpath_restraint         (int i_c, struct mtd_data_s *mtd_data);
 void pathref_findiff         (int i_c, struct mtd_data_s *mtd_data);
 void position_restraint      (int i_c, struct mtd_data_s *mtd_data);
 void elstpot_restraint       (int i_c, struct mtd_data_s *mtd_data);
 void puckering_restraint     (int i_c, struct mtd_data_s *mtd_data);
 void energy_restraint        (int i_c, struct mtd_data_s *mtd_data);
 void helix_restraint         (int i_c, struct mtd_data_s *mtd_data);
 void cmap_restraint          (int i_c, struct mtd_data_s *mtd_data);
 void alpharmsd_restraint     (int i_c, struct mtd_data_s *mtd_data);
 void antibetarmsd_restraint  (int i_c, struct mtd_data_s *mtd_data);
 void parabetarmsd_restraint  (int i_c, struct mtd_data_s *mtd_data);
 void sbernd_restraint        (int i_c, struct mtd_data_s *mtd_data);
 void zbernd_restraint        (int i_c, struct mtd_data_s *mtd_data);
 void pca_restraint           (int i_c, struct mtd_data_s *mtd_data);
 void poly_restraint          (int i_c, struct mtd_data_s *mtd_data);
 void func_restraint          (int i_c, struct mtd_data_s *mtd_data);
 void msd_restraint           (int i_c, struct mtd_data_s *mtd_data);
 void sprint_restraint        (int i_c, struct mtd_data_s *mtd_data);

 // bernd pathway  
 void init_bernd_evolution    ( struct sz_data **pmy_sz , int i_c);
 void do_bernd_evolution      ( struct sz_data *pmy_sz );
 void reparam_bernd_path      ( struct sz_data *pmy_sz );
 void dump_frameset_formatted  ( struct sz_data *pmy_sz );
 void dump_frameset_unformatted( struct sz_data *pmy_sz );

 // v1=v2+dr(v3-v4)
 void do_step_bernd           (  struct hybrid_frameset *v1,struct hybrid_frameset *v2, real **mat, real dr, struct hybrid_frameset *v3 ,  struct hybrid_frameset *v4 );
 #if defined (PLUMED_GROMACS) && defined (__cplusplus)
 extern "C" {
 #endif
 int  read_camshift      (char **word,int count,t_plumed_input *input,           FILE *fplog);
 void camshift_restraint      (int i_c, struct mtd_data_s *mtd_data);
 void camshiftens_restraint   (int i_c, struct mtd_data_s *mtd_data);
 #if defined (PLUMED_GROMACS) && defined (__cplusplus)
 } 
 #endif
 real coscut(real r,real r_0, real cut);
 real dcoscut(real r,real r_0, real cut);
 real tprod (rvec a,rvec b, rvec c );
 real generate_R(int i_c,struct mtd_data_s *mtd_data,rvec* R,rvec Rp,rvec Rdp,rvec RpxRdp);
 real puckering_Zeta (rvec R,rvec RpxRdp,real mod);
 void puckering_gradZeta(rvec gradZ,int index_i,int index_j,real* z,rvec* R,rvec Rp, rvec Rdp, rvec RpxRdp,real mod);
 real puckering_Q(real *z);
 void puckering_gradQ(rvec gradQ, real* z, real Q,rvec*R, rvec Rp, rvec Rdp, rvec RpxRdp,real mod,int index);
 real puckering_phi(real *z);
 void puckering_gradphi(rvec gradphi, real* z,rvec*R, rvec Rp, rvec Rdp, rvec RpxRdp,real mod,int index);
 real puckering_theta(real * z);
 void puckering_gradtheta(rvec gradtheta, real* z,rvec*R,rvec Rp,rvec Rdp,rvec RpxRdp,real mod,int index);
// GRID stuff
real spline(int ndim,real *dx,real *where,real *tabf,real *tabder,int* stride,real *der);
void grid_initialize(struct grid_s *grid);
void grid_resize_minigrid(struct grid_s *grid, real* delta, real cutoff);
void grid_create_one2multi(int **one2multi, int size, int ncv, int *bin);
void grid_addhills(struct grid_s *grid, real ww, real* ss,real* delta,int rank,int npe);
real grid_getstuff(struct grid_s *grid, real* ss0,real* force);
int  grid_multi2one(struct grid_s *grid, int* index_nd);
void grid_write_tofile(struct grid_s *grid);
void grid_read_fromfile(struct grid_s *grid, int bias); 
void grid_clone(struct grid_s *grid1, struct grid_s *grid2);
// misc routines
 void cmap_running(int i_c, struct cmap_inpack *inpack, struct cmap_pack *my_cmap_pack);
 void cmdist_eval(int i_c, int frame,struct cmap_inpack *inpack,struct cmap_outpack *outpack, 
      struct cmap_pack *my_cmap_pack,int dr1_calc);
 real pow2(real x);
 void power(real x,int p,int q,real *xp,real *xq);
 void read_sz_map(struct sz_data *my_sz, char file_maps[129], char file_maps2[129],
                        char file_group[129], int read_mapfile, FILE *fplog);
 int  read_sz_rmsd(struct sz_data *my_sz, FILE *fplog);
 int  read_sz_hybrid(struct sz_data *my_sz, FILE *fplog);
 int  read_sz_coord (char *filename, struct coordinates_frameset *p, FILE *fplog);
// coordination 
 void coord_newlist(int i_c, struct mtd_data_s *mtd_data);
 void coord_checklist(int i_c, struct mtd_data_s *mtd_data);
 void coord_restraint_nlist(int i_c, struct mtd_data_s *mtd_data);
 void coord_restraint_no_nlist(int i_c, struct mtd_data_s *mtd_data);
 // hbd copy/read simple
 int  hbd_read_simple (FILE *myfile, FILE *fplog, const char *string, struct hybrid_elem *, int hasfile);	
 void  hbd_copy_simple ( struct hybrid_elem *elem ); 
	void  hbd_dump_simple( struct hybrid_elem *elem , FILE *fp);
	void  hbd_dump_msd( struct hybrid_elem *elem , FILE *fp);

 //real hbd_vecmvec(struct hybrid_frameset *v1,  struct hybrid_frameset *v2, struct hybrid_frameset *v3, struct hybrid_frameset *v4 ,real **mat , real **dv1dcv,  real **dv2dcv,  real **dv3dcv,  real **dv4dcv, int allo);
 real hbd_vecmvec_ref(struct hybrid_frameset *v1,  struct hybrid_frameset *v2, 
								struct hybrid_frameset *v3,  struct hybrid_frameset *v4 ,
								struct hybrid_frameset *v5,  struct hybrid_frameset *v6 ,							
								real **mat , real *dv1dcv,  real *dv2dcv,  real *dv3dcv, 
								real *dv4dcv, real *dv5dcv,  real *dv6dcv, FILE *fplog);
 real hbd_distanddiff_ref(struct hybrid_frameset *v1,  struct hybrid_frameset *v2, 
								struct hybrid_frameset *v3,  						
								real **mat , real *dv1dcv,  real *dv2dcv,  real *dv3dcv,struct hybrid_frameset *d1,  
								  FILE *fplog);
void calc_diff_twoframes(struct hybrid_frameset *v1,  struct hybrid_frameset *v2, struct hybrid_frameset *v3,struct hybrid_frameset *d1 );
void check_hbd_vecmvec_ref(struct hybrid_frameset *r1,struct hybrid_frameset *r2,struct hybrid_frameset *r3,struct hybrid_frameset *r4,struct hybrid_frameset *r5,struct hybrid_frameset *r6,real **matrix, FILE  *fplog);
void clone_hybrid_frameset(struct hybrid_frameset **sink,struct hybrid_frameset *source,  int need_alloc_diff, FILE *fplog);
void destroy_hybrid_frameset(struct hybrid_frameset *sink, FILE *fplog);
void clone_hybrid_elem(struct hybrid_elem **sink,struct hybrid_elem *source,  int need_alloc_diff, FILE *fplog);
void destroy_hybrid_elem(struct hybrid_elem *sink, FILE *fplog);
void simple_elem_diff(struct hybrid_elem *diff,struct hybrid_elem  *e1,struct hybrid_elem  *e2,struct hybrid_elem  *e3,  FILE *fplog);
void msd_elem_diff(struct hybrid_elem *diff,struct hybrid_elem  *e1,struct hybrid_elem  *e2,struct hybrid_elem  *e3, FILE *fplog);
void check_msd_elem_diff(struct hybrid_elem *diff,struct hybrid_elem *v1, struct hybrid_elem *v2, struct hybrid_elem *v3, FILE *fplog);
//hbd msd
int   hbd_read_msd (  FILE *myfile, FILE *fplog, struct hybrid_elem *elem , int hasfile);
void  hbd_copy_msd ( struct hybrid_elem *elem  ) ;
// 
int read_pdb (struct pdb **mypdb, FILE *myfile, FILE *fplog);
void copy_pdb (struct pdb *pdb1, struct pdb *pdb2, FILE *fplog);
int allocate_pdb(struct pdb **mypdb,int i);
int deallocate_pdb(struct pdb *mypdb);
void readapt_rmsd_work_array(int i, FILE *fplog);
void readapt_rmsd_work_array_secondder(int i, FILE *fplog);
void allocate_rmsd_container(struct rmsd_container_t *rmsd_container,int i);
void copy_pdb_into_rmsd_container(struct pdb, struct rmsd_container_t *rmsd_container, FILE *fplog);
void msd_calculation_dynamic(struct rmsd_workstruct_s *work,int der_frameref_on, int do_rot, int do_center);
void msd_core_dynamic_norot(struct rmsd_workstruct_s *work, int do_center);
void msd_core_dynamic_simple(struct rmsd_workstruct_s *work,int do_center, int do_frameref_der, int do_rotmat_der);
void msd_core_dynamic_weighted(struct rmsd_workstruct_s *work,int do_center, int do_frameref_der);
void rmsd_dynamic_findiff_interface(struct rmsd_workstruct_s *work);
void clean_rmsd_work_array(struct rmsd_workstruct_s *work);
void calc_projector_test(struct sz_data *pmy_sz);	
void calc_intraframe_dist( struct sz_data *pmy_sz );
void calc_intraframe_diff( struct sz_data *pmy_sz );
void reparam_with_multiple_matrix( struct sz_data *pmy_sz , char *filename);
void calc_twoframe_dist( struct sz_data *pmy_sz , int first, int second,  char *matrixfile, FILE *fplog);
void point_to_matrix(int i,struct sz_data *my_sz);
// int  hbd_copy_target ( struct hybrid_elem *elem ); 
// real hbd_metrics_target ( struct hybrid_elem *run,  struct hybrid_elem *ref );
 int  hbd_collect_config ( struct hybrid_frameset *running  );
 int  hbd_collect_jacobian ( struct hybrid_frameset *running, real ** mathybrid ,  real ** myinverse, FILE *fplog , int absval);
// int  hbd_metrics ( struct hybrid_frameset *running , struct hybrid_frameset *reference , struct cmap_outpack *outpack, real **mat); 
 int  hbd_metrics_new ( struct hybrid_frameset *running , struct hybrid_frameset *reference , struct cmap_outpack *outpack, real **mat ,FILE *fplog); 
 void  test_hbd_metrics_new ( struct hybrid_frameset *running , struct hybrid_frameset *reference , struct cmap_outpack *outpack, real **matrix, FILE  *fplog);
 void  spath_restraint_testder(int i_c, struct mtd_data_s *mtd_data);
 void  zpath_restraint_testder(int i_c, struct mtd_data_s *mtd_data);
 void msd_calculation(struct coordinates_frameset *pframeset,struct cmap_inpack *c_inpack,
                             struct cmap_outpack *c_outpack,real dmsd_dr1[3][MAXATOMS_PATH],int der_frameref_on, int norot, int nocenter);
 int  rmsd_pack(struct rmsd_inpack inpack,struct rmsd_outpack *outpack,int iopt,int iopt2);
 int  rmsd_mini_pack(struct rmsd_inpack inpack,struct rmsd_mini_outpack *outpack,int iopt,int simple, int permissive);
 int  rmsd_mini_pack_fake(struct rmsd_inpack inpack,struct rmsd_mini_outpack *outpack, int nocenter, int simple);
 int  rmsd_findiff_interface(struct rmsd_inpack inpack,struct rmsd_mini_outpack *outpack);
 void mean_rmsd(struct sz_data *pmy_sz, real dCV_dr1[MAXFRAMES_PATH][3][MAXATOMS_PATH],
                        int i_c, FILE *fplog); 
 void mean_map(struct sz_data *pmy_sz, real dCV_dcm[MAXFRAMES_PATH][MAXDIM_CMAP],
                        int i_c, FILE *fplog);
 void dmsd_calculation(int i_c,struct coordinates_frameset *pframeset,struct cmap_inpack *c_inpack,
                             struct cmap_outpack *c_outpack,real dmsd_dr1[3][MAXATOMS_PATH]);
 // stuff for gradient projection
 void couple2list( struct coupling_ll **atlist ,int *at1,int nat1,int *at2,int nat2);
 void freecouple2list( struct coupling_ll **atlist);
 void scancouple( struct coupling_ll *atlist);
 void setup_projections(struct proj_grad_s * );   
 void calc_projections(struct proj_grad_s *);   
// misc
 void cite_please (const char* re, FILE *fplog);
 void disclaimer (FILE *fplog);
 // allocators and destructors
 real *float_1d_array_alloc(int ii);
 real **float_2d_array_alloc(int ii,int jj);
 real ***float_3d_array_alloc(int i,int j,int k);
 real ****float_4d_array_alloc(int i,int j,int k,int l);
 int  *int_1d_array_alloc(int ii);
 int  **int_2d_array_alloc(int ii,int jj);
 int free_1dr_array_alloc(real *xx); // 1d real
 int free_2dr_array_alloc(real **xx,int ii); // 2d real
 int free_3dr_array_alloc(real ***xx,int ii,int jj); // 3d read
 int free_4dr_array_alloc(real ****xx,int ii,int jj, int kk); //4d real
 int free_1di_array_alloc(int *xx);  //1d integer
 int free_2di_array_alloc(int **xx,int ii);  //2d integer
 // neighbour list tools for quicksorting
 void realquicksort ( real *v , int *ind , int left , int right ); // quicksort for neighbour list search
 void swap (real *v,int *ind,int i,int j); // used by quicksort
 int ql77_driver(real m[][4],real* lambda);
 void rank3_to_ql77(double in[3][3], double evector[3][3], double evalue[3]);
 int ql77 (int n,double *x,double *d);	
 void import_ow(real **weight ,char* filename , int start, int end , int *indexes);
 int matinverse (real **a, real **id, int n, FILE *fplog);
 void parse_fixwidth(char *str,int i, int j, char *str2);
//! RECONMETA INCLUDES  Gareth Tribello
#ifdef RECONMETAD 
//   #include "recon_cbind.h"
   MYEXT int reconOn;
   MYEXT struct recon_data_s  reconinpt;
   MYEXT void* myreconObj;
#endif
#ifdef CVS
   MYEXT struct bespoke_data_s bespoke_input;
   MYEXT void* mybespokeObj;
#endif

/* COMMON DATA STRUCTURES  */
    MYEXT struct mtd_data_s    mtd_data;
    MYEXT struct dafed_s       dafed[nconst_max];		// d-AFED #####
    MYEXT struct dafed_control_s       dafed_control; // d-AFED #####
    MYEXT struct steer         cvsteer;
    MYEXT struct tamd          tamd;
    MYEXT struct abmd          abmd;
    MYEXT struct constraint    cvcnstr;
    MYEXT struct wall          cvw;
    MYEXT struct logical_s     logical;
    MYEXT struct colvar_s      colvar;
    MYEXT struct hills_s       hills;
    MYEXT struct steerplan_s   steerplan;
    MYEXT struct nlist_s       nlist[nconst_max];
    MYEXT int   firstTime;					// first PluMed step
    MYEXT real   Vrecon;                                        // Reconaissance metadynamics potential
    MYEXT real   Vext;                                          // External potential
    MYEXT real   Vconstr;                                       // Constraint potential
    MYEXT real   Vsteerplan;                                    // Steerplan potential
    MYEXT real   fext[nconst_max];                              // External forces
// path stuff
    MYEXT struct sz_data my_sz_list[NMAX_PATH];
	MYEXT struct rmsd_workstruct_s rmsd_workstruct;  // a general rmsd container that is allocated/reallocated on the fly
	MYEXT struct couplingmatrix_s couplingmatrix;
    MYEXT int ic_to_sz[nconst_max];
    MYEXT int nsz;
    MYEXT int kill_me[nconst_max];
// GRID stuff
    MYEXT struct grid_s       bias_grid;
// COMMIT
    MYEXT struct commit_s   commit;
// alpha/betarmsd stuff
    MYEXT struct ref_dist_mat_s  ref_dist_mat;
// sprint stuff
    MYEXT struct sprint_data_s sprint_data;
// external potential
    MYEXT struct grid_s     extpot;
//target distribution, if being used
    MYEXT struct grid_s    target_grid;
// camshift stuff
    MYEXT struct cam_shift_s cam_shift;
// stopwhen stuff
    MYEXT struct stopwhen_s stopwhen;
// interval staff
    MYEXT struct interval_s cvint;
// function stuff
    MYEXT struct  mathfunction_s      mathfunction[nconst_max];

#if defined(DL_POLY) || defined(AMBER) || defined(PLUMED_QESPRESSO)
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

#if defined (NAMD) || defined (LAMMPS_PLUMED) || defined (PLUMED_CPMD)
inline real min(real a, real b) { if(a < b){ return a;}else{b;}; };
#endif

#ifdef NAMD
};
#elif LAMMPS_PLUMED
};
#elif PLUMED_CPMD
// we need this because Plumed is a class
};
// this makes the Plumed object global
extern Plumed Plumed_simulation;
#endif

