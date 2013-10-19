#ifdef RECONMETAD  
#include "recon_basins.h"

namespace rp {

std::istream& operator>>(std::istream& istr, basinObj& bas){
   istr>>bas.refval>>bas.covar>>bas.diffConst1;                // Read in the line
   return istr;
}

std::ostream& operator<<(std::ostream& ostr,const basinObj& bas){
   ostr<<" "<<bas.refval<<" "<<bas.covar<<" "<<bas.diffConst1;     
   return ostr; 
}  

void basinObj::setBasin( const index& periodic, const rVector& center, const rMatrix& cov, const real& sizeParam, const real& w, const real& h ){
   refval=center;                                     // Set the positiion of the center
   covar=cov;                                         // The covariance
   invertMat( covar, metric );                        // The metric (inverse of covariance) 
   det=logdet( covar );                               // The logarithm of the determinant of the covariance
   initPsize=Psize=sqrt( ncv-1.0 ) + sizeParam;       // The initial size of the basin 
   maxSize=0.0;                                       // Set the maximum for expansion of periodic variables
   if(periodic==1){ for(index i=0;i<ncv;++i) maxSize+=metric(i,i); } maxSize=sqrt(maxSize);
//   std::pout<<"Setting maximum size of this basin equal to "<<maxSize<<std::endl;
   nonions=0;                                         // The number of onions
   width=w;                                           // The width of the onions
   height=h;                                          // The height of the onions
}

void basinObj::restartBasin( const index& periodic, const real& sizeParam, const real& w, const real& h ){
   invertMat( covar, metric );                      // The metric (inverse of covariance)
   det=logdet( covar );                             // The logarithm of the determinant of the covariance
   initPsize=Psize=sqrt( ncv-1.0 ) + sizeParam;     // The initial size of the basin 
   maxSize=0.0;                                     // Set the maximum for expansion of periodic variables
   if(periodic==1){ for(index i=0;i<ncv;++i) maxSize+=metric(i,i); } maxSize=sqrt(maxSize);
   nonions=0;                                       // The number of onions
   width=w;                                         // The width of the onions
   height=h;                                        // The height of the onions
//   diffConst2=expandParam - (diffConst1*econst) / ( initPsize );  // The second parameter for basin expansion
}

// These caclulate the (square of the) distance (in the space of the metric) between two arbitrary points.
real basinObj::calcD( const rVector& cv1, const rVector cv2 ) const {
   rVector tmp( cv1 ), derivs( ncv ); real d; tmp-=cv2;
   matVecMult( metric, tmp, derivs ); d=dotProduct( tmp, derivs ) ;
   return d;
}

real basinObj::calcD( const rVector& period, const rVector& cv1, const rVector cv2 ) const {
   // Find the difference from refval
   rVector tmp( cv1 ); tmp-=cv2;
   // Compute sines and cosines of the vectors
   rVector costmp( ncv ), sintmp( ncv ), sinout( ncv );
   for(index i=0;i<ncv;++i){
      tmp[i]*=period[i]; sintmp[i]=sin( tmp[i] ); costmp[i]=cos( tmp[i] );
   }

   // Now calculate d 
   real d=0.0;
   for(index i=0;i<ncv;++i){ d+=2*( 1 - costmp[i] )*metric(i,i); }

   d+=dotProduct( sintmp, sinout);
   if( d<0.0 ){ 
      WARNING("THE DISTANCE FROM THE METRIC IS LESS THAN ZERO SUGGESTING AN ILL CONDITIONED METRIC");
      d=fabs(d);
   }

   return d;
}

// These calculate the distance (in the space of the metric between a point a the center of the basin
real basinObj::calcR( const rVector& colvars, rVector& derivs ) const {
   rVector tmp( colvars ); tmp-=refval;
   real d; matVecMult( metric, tmp, derivs ); d=sqrt( dotProduct( tmp, derivs ) );
   derivs/=d;
   return d;
}

real basinObj::calcR( const rVector& period, const rVector& colvars, rVector& derivs ) const {

   // Find the difference from refval
   rVector tmp( colvars ); tmp-=refval;
   // Compute sines and cosines of the vectors
   rVector costmp( ncv ), sintmp( ncv ), sinout( ncv );
   for(index i=0;i<ncv;++i){
      tmp[i]*=period[i]; sintmp[i]=sin( tmp[i] ); costmp[i]=cos( tmp[i] );
   }
 
   // Now calculate R and derivatives
   real d=0.0; sinout=0.0; 
   for(index i=0;i<ncv;++i){ 
     d+=2*( 1 - costmp[i] )*metric(i,i);
     for(index j=0;j<ncv;++j){
        if(i!=j){ sinout[i]+=metric(i,j)*sintmp[j]; }
     }
     derivs[i]=period[i]*( metric(i,i)*sintmp[i] + sinout[i]*costmp[i] );
   }

   d+=dotProduct( sintmp, sinout);
   if( d<0.0 ){ 
      WARNING("THE DISTANCE FROM THE METRIC IS LESS THAN ZERO SUGGESTING AN ILL CONDITIONED METRIC");
      d=fabs(d);
   }
   d=sqrt(d); derivs/=d;

   return d;
}

real basinObj::calcPotential( const index& ncv, const rVector& colvars, const rVector& period, rVector& derivatives ){

  // Is the data periodic? Calculate position with respect to center and the derivatives in manner accordingly
  rVector derivs( ncv );
  if(period[0]!=0){ R=calcR( period, colvars, derivs);  }
  else{ R=calcR( colvars, derivs ); }

  // Don't compute the onions if R is very large
  if( R>BASIN_CUTOFF*Psize ){ return 0.; }
  
  // Now sum the contributions to the energies and forces from all the onions
  real eng=0.0, F=0.0; 
  // This calculates the potential due to the onions
  for(index i=0;i<nonions;++i){ eng+=onionList[i].feelOnion( R, F ); }

  // Update the derivatives
  for(index i=0;i<ncv;++i){ derivatives[i]+=F*derivs[i]; }

  return eng;
}

void basinObj::addOnion( mpiostream *onionfile, const index& bas, const real& time ){
  onionObj newOnion( R, width, height );
  onionList.push_back( newOnion ); nonions++; 
  (*onionfile)<<time<<" "<<bas<<" "<<newOnion<<" "<<Psize<<std::endl;
  return;
//  return newOnion.heightAtZero();
}

void basinObj::inflateBasin( mpiostream *onionfile, const real& econst, const index& bas, const real& time ){

  // This controls basin expansion
  real Pexpand; Pexpand = (econst*diffConst1/Psize);                // + diffConst2;
  real rand; rand=drand48();
  if( rand<Pexpand ){
     // This prevents overexpansion in periodic variables 
     if( maxSize==0.0 || (Psize+width)<maxSize ){   
        Psize+=width;
        onionObj newOnion( 0.0,0.0,0.0 );
        (*onionfile)<<time<<" "<<bas<<" "<<newOnion<<" "<<Psize<<std::endl;
     }
  } 
}

std::ostream& operator<<(std::ostream& ostr,const onionObj& o ){
   ostr<<o.position<<" "<<o.width<<" "<<o.height;
   return ostr;
}

std::istream& operator>>(std::istream& istr, onionObj& o ){
   istr>>o.position>>o.width>>o.height;
   return istr;
}

real onionObj::feelOnion( const real& r, real& fp ) const {
  real pot, tmp1, tmp2, dv1, dv2;
  dv1=r-position; tmp1=-0.5*pow(dv1/width,2); 
  dv2=r+position; tmp2=-0.5*pow(dv2/width,2);
  
  // Calculate the potential due to this onion
  tmp1=height*exp(tmp1); tmp2=height*exp(tmp2); pot=tmp1+tmp2;

  // Calculate the force due to this onion
  fp+=tmp1*( dv1/pow(width,2) ) + tmp2*( dv2/pow(width,2) );

  return pot;
}

index closerBasin( const basinObj& bas1, const basinObj& bas2 ){
   if( bas1.R<bas2.R ) return 1;
   else if( bas2.R<bas1.R ) return 2;
   return 1;
}

real overlap( const index& ncv, const rVector& period, const basinObj& clus , const basinObj& bas ){
   rMatrix fullCovar(ncv, ncv), fullMetric(ncv,ncv);

   // Get the full covariance
   fullCovar=bas.covar; fullCovar*=bas.Psize/bas.initPsize; 
   fullMetric=clus.covar; fullMetric*=clus.Psize/clus.initPsize; 
   fullCovar+=fullMetric;
   // Invert the covariance
   invertMat( fullCovar, fullMetric );

   // Calculate the vector connecting the two basins
   rVector tmp( clus.refval ); tmp-=bas.refval; pbc( tmp, period ); 

   // Now calculate the exponent
   real exponent=0.0;
   for(index i=0;i<ncv;++i){
       for(index j=0;j<ncv;++j){
          exponent+=tmp[i]*fullMetric(i,j)*tmp[j]; 
       }
   }

   real overl; overl=exp(-(exponent/4) );

   // Ensure that the determinant is corrected with the size
   real det1; det1=log((bas.Psize/bas.initPsize))*ncv +bas.det;
   real det2; det2=log((clus.Psize/clus.initPsize))*ncv +clus.det;
   real det12; det12=logdet( fullCovar );

   // And finally calculate Q - using logs to prevent underflow
   double Q; Q = log(2) * double(ncv) / 2.0 + det2/4. +det1/4. - det12/2.; Q=exp(Q);

#ifdef DEBUG
   if( Q*overl>1.05 || Q*overl<0.0 ) ERROR("Dubious overlap between basins "<<Q*overl);
#endif

   return Q*overl;
}

// Calculate the distace between two basins as the arithmetric average of the distances in the metric
real basinDistance( const index& ncv, const rVector& period, const basinObj& clus , const basinObj& bas ){
   real dist=0;

   // Convert period to periodicity
   rp::rVector periodicities(ncv); periodicities=0;
   if( period[0]!=0) {
     for(rp::index i=0;i<ncv;++i) periodicities[i]= 2*R_PI/period[i];
   }

   // Get distance of clus from bas 
   rp::rVector colvars( ncv ), derivs( ncv );
   clus.getCenter( colvars ); 
   if( periodicities[0]!=0 ){ dist=bas.calcR( periodicities, colvars, derivs ); }
   else{ dist=bas.calcR( colvars, derivs ); }

   // Get distance of bas from clus
   bas.getCenter( colvars ); 
   if( periodicities[0]!=0 ){ dist+=clus.calcR( periodicities, colvars, derivs ); }
   else{ dist+=clus.calcR( colvars, derivs ); }

   return dist/2.0;
}

//#ifdef RECON_DRIVER
//// Test if we are in basin
//index basinObj::inMetadRegion( const rVector& period, const rVector& colvars ) const {
//   rVector derivs( ncv );
//   if(period[0]!=0){ if( calcR( period, colvars, derivs)<Psize){return 1; } }
//   else if( calcR( colvars, derivs)<Psize ){ return 1; }
//   return 0;
//}
//// Distance between two basins
//real dist( const index& ncv, const rVector& period, const basinObj& bas1, const basinObj& bas2 ){
//#ifdef DEBUG
//   if( bas1.refval.size()!=bas2.refval.size() ) ERROR("SIZE MISMATCH");
//#endif
//   // Establish difference between vectors
//   rVector diff( bas1.refval ); diff-=bas2.refval; pbc( diff, period );
//   
//   return sqrt( dotProduct(diff, diff) );
//}
//
//real basinObj::dist( const rVector& period, const rVector& colvars ) const {
//#ifdef DEBUG
//    if( colvars.size()!=refval.size() ) ERROR("SIZE MISMATCH");
//#endif
//
//    // Establish difference between vectors
//    rVector diff( colvars ); diff-=refval; pbc( diff,period );
//
//    return sqrt( dotProduct(diff, diff) );
//}
//
//real basinObj::calcBias( const real& d ) const {
//    // Now sum the contributions to the energies and forces from all the onions
//    real eng=0.0, F=0.0;
//    // This calculates the potential due to the onions
//    for(index i=0;i<nonions;++i){ eng+=onionList[i].feelOnion( d, F ); }
//    return eng;
//}
//#endif

// These are the restart functions for reading in basins and onions
void readBasins( const index& ncv, const index& periodic, const real& sizeParam, const real& width, const real& height, std::vector<basinObj>& basinList ){
    // Open the basin file
    std::ifstream* ibasinfile=new std::ifstream("BASINS");
    if (!ibasinfile->good()) ERROR("Something bad happened while opening BASINS file for restart");

    basinObj newbasin(ncv); real tstep, extranum; index num; std::string line;
    while ( getline( (*ibasinfile), line) ) {
       // Convert the read in line to a stringrstream and read in
       std::stringstream input (line);   // input >> num >> tstep >> newbasin; 
       if( ( input >> num >> tstep >> newbasin ).fail() ) ERROR("WRONG FORMAT IN BASIN FILE");
       // This checks there is nothing else
       if( !(input>>extranum).fail() ) WARNING("STRANGE FORMAT IN BASIN FILE"); 
       // Restart the basin
       newbasin.restartBasin( periodic, sizeParam, width, height );
       // Add it to the basin list
       basinList.push_back(newbasin);
    }
    ibasinfile->close();
//    std::pout<<"#Read "<<basinList.size()<<" basins during restart"<<std::endl;
}

index readOnions( std::vector<basinObj>& basinList ){
    // Open the file
    std::ifstream* ionionfile=new std::ifstream("ONIONS");
    if (!ionionfile->good()) ERROR("Something bad happened while opening ONIONS file for restart");

    index bas,n=0; real tstep; std::string line; 
    onionObj newOnion; real size, extranum; // real conserved=0.;
    while( getline( (*ionionfile), line) ) {
       // Convert the line into a stringstream and read in
       std::stringstream input (line);  // input >> tstep >> bas >> newOnion >> size;
       if( (input >> tstep >> bas >> newOnion >> size).fail() ) ERROR("WRONG FORMAT IN ONION FILE");
       // Check on format
       if ( !(input>>extranum).fail() ) WARNING("STRANGE FORMAT FOR ONION FILE");
       // Add Onion to list
       basinList[bas].restartOnion( newOnion, size ); n++; 
       //conserved+=basinList[bas].restartOnion( newOnion, size ); n++;
    }
//    std::pout<<"#Read in "<<n<<" onions during restart"<<std::endl;
    ionionfile->close();
    return n;
}

}  // END OF NAMESPACE
#endif
