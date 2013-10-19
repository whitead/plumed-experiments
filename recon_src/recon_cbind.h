#ifndef _RECON_CBIND
#define _RECON_CBIND 1

#include <stdio.h>

#define nscales_max 10

struct recon_data_s{
 // Are we doing basin monitoring only
 // int monitor;
 // Data on the cvs we are using for reconnaissance
 int nconst;
 int* cvlist;
 // Data for clustering
 int nscales, fuzzy;
 int runFreq[nscales_max]; 
 int storeFreq[nscales_max];
 // Data for hills
 int stride;
 double height;
 double width;
// double temp; 
 // Data for controlling basins
 double basTol;
 double diffconst;
 double iSize; 
 // Logical stuff to control restart
 int restart;
};

#ifdef CVS
#define MAX_BESPOKE 2
#define BESPOKE_MAXFILENAME 1024
struct bespoke_data_s{
 int nbespoke;
 double smooth;
 int nspline[MAX_BESPOKE];
 int ngrid[MAX_BESPOKE];
 char filename[BESPOKE_MAXFILENAME];
};
#endif

#ifdef __cplusplus 
extern "C" {
#endif

#ifdef CVS
// These are for the interface with bespoke collective coordinates
void create_bespoke( int ncv, int nbespoke, void **);
void setup_bespoke( struct bespoke_data_s inp, int ncv, int nbespoke, double* periods, void *bespoke, FILE *fplog );
void calculate_bespoke_cvs( int ncv, int dimout, double* cv_in, double* cv_out, double* derivatives, void *bespoke);
//double calculate_bespoke_error( int ncv, double* cv_in, void *bespoke);
double calculate_bespoke_error(void *bespoke);
#endif

// These are the interface with reconnaissance metadynamics
#ifndef DRIVER
void create_recon(void **);
void setup_recon(double* periods, double tstep, struct recon_data_s reconinpt, void *recon, FILE *fplog);
double dometa_recon(int timestep, int ncv, double* ss0, double* recon_der, void *recon);
// int dometa_monitor(int ncv, double* ss0, void *recon);
// double getCons_recon( void* recon );

// These are the routines for testing the derivatives
void setup_test_recon( int ncv, double* ss0, void* recon );
void recon_calc_cv( int ncv, double* ss0, double* pos, void* recon);
void recon_derivatives( int ncv, double* ss0, double* der, void* recon );
#endif

#ifdef __cplusplus 
}
#endif


#endif //ends #ifndef _RECON_CBIND
