#ifdef RECONMETAD 

#ifndef DRIVER
#include "recon_metad.h"
#endif

#ifdef CVS
#include "recon_bespoke.h"
#endif

extern "C" {
#include "recon_cbind.h"

// Bindings for bespoke collective coordinates

#ifdef CVS
void create_bespoke( int ncv, int nbespoke, void **bespoke ){
   if( nbespoke!=MAX_BESPOKE ) ERROR("WRONG NUMBER OF BESPOKE CVS IN INPUT - SHOULD BE "<<MAX_BESPOKE);
   rp::index ncolvars, dimout; 
   ncolvars=ncv; dimout=nbespoke; 
   *bespoke= new rp::bespokeCVObj( ncolvars,dimout );
}

void setup_bespoke( struct bespoke_data_s inp, int ncv, int nbespoke, double* periods, void *bespoke, FILE *fplog ){ 
   // Get the periods
   rp::rVector cv_per(ncv); for(int i=0;i<ncv;++i) cv_per[i]=periods[i]; 
   std::string ifile; ifile=inp.filename;   // Get the file
   rp::real smooth; smooth=inp.smooth;      // Get the smoothing parameter
   // Get the grid for spline
   std::vector<rp::index> nspline(MAX_BESPOKE); for(int i=0;i<nbespoke;++i) nspline[i]=inp.nspline[i]; 
   // Get the grid for integration
   std::vector<rp::index> ngrid(MAX_BESPOKE); for(int i=0;i<nbespoke;++i) ngrid[i]=inp.ngrid[i];
   // Call the setup routine
   ((rp::bespokeCVObj*)bespoke)->setup( ifile,cv_per,smooth,nspline,ngrid,fplog ); 
}

void calculate_bespoke_cvs( int ncv, int dimout, double* cv_in, double* cv_out, double* derivatives, void *bespoke){
   rp::rVector colvars(ncv), embededP(dimout); rp::rMatrix der( dimout,ncv );
   for(int i=0;i<ncv;++i){ colvars[i]=cv_in[i]; }
   ((rp::bespokeCVObj*)bespoke)->embedPoint( colvars, embededP, der );
   int k=0;
   for(int i=0;i<dimout;++i){
      cv_out[i]=embededP[i];
      for(int j=0;j<ncv;++j){ derivatives[k]=der(i,j); k++; }
   }
}

//double calculate_bespoke_error( int ncv, double* cv_in, void *bespoke){
double calculate_bespoke_error( void *bespoke) {
//   rp::rVector colvars(ncv); double error;
//   for(int i=0;i<ncv;++i){ colvars[i]=cv_in[i]; }
//   error=((rp::bespokeCVObj*)bespoke)->reconstructError( colvars );
   double error; error=((rp::bespokeCVObj*)bespoke)->reconstructError();
   return error;
}
#endif

// Bindings for reconnaissance metadynamics 

#ifndef DRIVER 
void create_recon(void ** recon ){  
   *recon= new rp::reconObj(); 
}
 
void setup_recon(double* periods, double tstep, struct recon_data_s reconinpt, void *recon, FILE *fplog){ 
   rp::rVector cv_per( reconinpt.nconst );
   for(int i=0;i<reconinpt.nconst;++i) cv_per[i]=periods[i];
   ((rp::reconObj*)recon)->setup(reconinpt.nconst,cv_per,tstep,reconinpt, fplog);

//   if( reconinpt.monitor==1 ){ ((rp::reconObj*)recon)->setupMonitor(reconinpt.nconst,cv_per); }
//   else{ ((rp::reconObj*)recon)->setup(reconinpt.nconst,cv_per,tstep,reconinpt, fplog); } 
}

double dometa_recon(int timestep, int ncv, double* ss0, double* recon_der, void *recon){

   rp::rVector colvars(ncv), derivatives(ncv); rp::real ene;
   for(int i=0;i<ncv;++i) colvars[i]=ss0[i];
   ene=((rp::reconObj*)recon)->doReconMetad(timestep, colvars, derivatives);
   for(int i=0;i<ncv;++i) recon_der[i]=derivatives[i];
   return ene;
}

//int dometa_monitor(int ncv, double* ss0, void *recon){
//   rp::rVector colvars(ncv); int current_basin;
//   for(int i=0;i<ncv;++i) colvars[i]=ss0[i];
//   current_basin=((rp::reconObj*)recon)->doMonitor( colvars );
//   return current_basin;
//} 

// double getCons_recon( void *recon ){
//   rp::real ene; ene=((rp::reconObj*)recon)->getConservedQuantity();
//   return ene;
// }
#endif

}
#endif
