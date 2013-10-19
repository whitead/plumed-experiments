#ifdef RECONMETAD 
#include "recon_metad.h"

extern "C" {
#include "recon_cbind.h"

void setup_test_recon( int ncv, double* ss0,  void* recon ){
  rp::rVector colvars(ncv);
  for(int i=0;i<ncv;++i) colvars[i]=ss0[i];
  ((rp::reconObj*)recon)->setup_test_derivatives( colvars );
}

void recon_calc_cv( int ncv, double* ss0, double* pos, void* recon){

   rp::rVector colvars(ncv), derivatives(ncv); rp::real r;
   for(int i=0;i<ncv;++i) colvars[i]=ss0[i];
   ((rp::reconObj*)recon)->test_derivatives_calc_cv(colvars, r );
   (*pos)=r;

}

void recon_derivatives( int ncv, double* ss0, double* der, void* recon ){

   rp::rVector colvars(ncv), derivatives; // , Qderivatives; 
   for(int i=0;i<ncv;++i) colvars[i]=ss0[i];
   ((rp::reconObj*)recon)->test_derivatives_calc_derivs(colvars, derivatives );
   for(int i=0;i<ncv;++i){ der[i]=derivatives[i]; } // Qder[i]=Qderivatives[i]; }
}

}

namespace rp {

// This routine creates a random basin
void reconObj::setup_test_derivatives( const rVector& colvars ){

   std::pout<<"Setting up test derivatives "<<ncv<<std::endl;

   rVector center( ncv ); rMatrix covar( ncv, ncv ), tmp( ncv, ncv);

   // Create a random center and random covariance
   covar=0.0;
   for(index i=0;i<ncv;++i){ 
     center[i]=colvars[i]+(drand48()-0.5)*0.1;
     for(index j=0;j<ncv;++j){tmp(i,j)=drand48(); }
   }
   tmatMatMult( tmp, tmp, covar );

   std::pout<<"Check center 1 "<<colvars<<std::endl;
   std::pout<<"Check center 2 "<<center<<std::endl;

//   // Now create a metric
//   index dim; dim=int( ncv / 2 );
//   metric thismetric( ncv ); 
//   // thismetric.setup( ncv, covar ); 
//   thismetric.setup( dim, covar );

   // Now create a basin and put it in the basin list
   basinObj mybasin(ncv );  
   index periodic_data=0; if(period[0]!=0.0){periodic_data=1;}
   mybasin.setBasin( periodic_data, center, covar, sizeParam, width, height );
   basinList.push_back(mybasin); 
   std::pout<<"basinList[0] "<<basinList[0]<<std::endl;
}

void reconObj::test_derivatives_calc_cv( const rVector& colvars, real& pos ){
   if(basinList.size()>1) ERROR("IN TEST DERIVATIVES THERE IS MORE THAN ONE BASIN "<<basinList.size());

   rVector derivs( ncv );
   if(period[0]==0.){pos=basinList[0].calcR( colvars, derivs );}
   else{pos=basinList[0].calcR( period, colvars, derivs );}
//   Ppos=basinList[0].calcP( colvars, shit ); 
//   Qpos=basinList[0].calcQ( colvars, shit ); 
}

void reconObj::test_derivatives_calc_derivs( const rVector& colvars, rVector& derivs ){
   if(basinList.size()>1) ERROR("IN TEST DERIVATIVES THERE IS MORE THAN ONE BASIN");

   real dum; derivs.resize(ncv);
   if(period[0]==0.){dum=basinList[0].calcR( colvars, derivs );}
   else{dum=basinList[0].calcR( period, colvars, derivs );}
//   shit=basinList[0].calcQ( colvars, Qderivs );
}

// real basinObj::calcQ( const rVector& colvars, rVector& derivs ) const {
//    rVector tmp( ncv ); 
//    tmp=colvars; tmp-=refval;
//    real R; R=mymetric.Qtransfer( tmp, derivs );
//    return R;
// }

//void reconObj::setup_test_onions( const rVector& colvars ){
//
//   std::pout<<"Setting up test onions "<<ncv<<std::endl;
//
//   rVector center( ncv ); rMatrix covar( ncv, ncv );
//
//   // Create a center and a random diagonal covariance
//   covar=0.0;
//   for(index i=0;i<ncv;++i){ center[i]=colvars[i]; covar(i,i)=drand48(); }
//
//   // Now create a metric
//   index dim; dim=int( ncv / 2 );
//   metric thismetric( ncv ); 
//   // thismetric.setup( dim, covar );
//   thismetric.setup( ncv, covar );
//
//   std::pout<<"Created a metric"<<std::endl;
//
//   // Now create a basin and put it in the basin list
//   basinObj mybasin(ncv);
//   mybasin.setBasin( center, covar, sizeParam );
//   // Only the height and width are important in setting this and the qsize (don't know how to set this)
//   mybasin.setParams( width, height );
//   basinList.push_back(mybasin);
//   std::pout<<"My basin "<<mybasin<<std::endl;
//  
//   // Create a few random onions
//   for(index i=0;i<10;++i){
//     real tmp,pos; tmp=sqrt( (ncv/2) - 1 ) + 3.0; 
//     pos=tmp*drand48();
//     std::pout<<"Created an onion at "<<pos<<std::endl;
//     basinList[0].addOnion( pos ); 
//   }
//   nbasins=basinList.size();
//}
//
//void reconObj::test_onions( rVector& colvars ){
//   real teststep = 0.001, invstep = 1000., highene, lowene, eng, ene, numforce;
//   rVector derivatives( ncv );
//
//   std::pout<<"In test "<<colvars<<std::endl;
//
//   ene=0.0;
//   for(index i=0;i<ncv;++i){
//      colvars[i]+=teststep; highene=0.0;
//      for(index j=0;j<nbasins;++j) highene+=basinList[j].calcPotential( ncv, ene, colvars, period, derivatives );
//      colvars[i]-=2*teststep; lowene=0.0;
//      for(index j=0;j<nbasins;++j) lowene+=basinList[j].calcPotential( ncv, ene, colvars, period, derivatives );
//      colvars[i]+=teststep; eng=0.0; derivatives=0.0;
//      for(index j=0;j<nbasins;++j) eng+=basinList[j].calcPotential( ncv, ene, colvars, period, derivatives );
//      numforce=0.5*((highene*invstep)-(lowene*invstep));
//      std::pout<<"FORCE ON CV "<<i<<" ANALYTICAL "<<derivatives[i]<<" NUMERICAL "<<-numforce<<" DIFF "<<derivatives[i]-numforce<<std::endl; 
//   }
//   abort();
//} 
//
//void basinObj::addOnion( const real& Rp ){
//   onionObj newOnion(Rp, width, height ); // Create a new onion
//   onionList.push_back( newOnion );       // Push it back to the onion list
//}

}
#endif
