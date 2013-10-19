#ifndef __rp_basins
#define __rp_basins

#include <vector>
#include <valarray>
#include <iostream>
#include <string>
#include <fstream>

#include "recon_types.h"
#include "recon_utils.h"
//#include "recon_metric.h"

#define BASIN_CUTOFF 5.0

namespace rp {

class onionObj {
friend std::istream& operator>>(std::istream&,onionObj&);
friend std::ostream& operator<<(std::ostream& ,const onionObj& );
private:
    real width,position,height;
public:
    onionObj() : width(0), position(0), height(0) {}
    onionObj( const real& pos, const real& w, const real& h ) : width(w), position(pos), height(h) {}
    onionObj( const onionObj& o ) : width(o.width), position(o.position), height(o.height) {}
    real feelOnion( const real& r, real& f ) const;
    index hasHeight() const {if(height>0.0){return 1;} return 0;}
//    real heightAtZero() const { return height*(1+exp(-(2.*position*position) / (width*width) ) ); }
};

std::istream& operator>>(std::istream& istr,onionObj& o);
std::ostream& operator<<(std::ostream& ostr,const onionObj& o );

class basinObj {
friend std::istream& operator>>(std::istream&,basinObj&);
friend std::ostream& operator<<(std::ostream&,const basinObj&);
friend real overlap( const index& ncv, const rVector&, const basinObj& , const basinObj& );
friend index closerBasin( const basinObj& , const basinObj& );
friend real basinDistance( const index& , const rVector& , const basinObj& , const basinObj& );
//#ifdef RECON_DRIVER
//friend real dist( const index& ncv, const rVector&, const basinObj& , const basinObj& );
//#endif
private:
    index ncv;
    index nonions;
    real det, R;
    rVector refval;
    rMatrix covar, metric;
    real Psize, initPsize, diffConst1, maxSize, width, height;
    std::vector<onionObj> onionList;
public:
    basinObj(index n=0) : ncv(n), nonions(0), det(0), R(0.), refval(ncv), covar(ncv,ncv), metric(ncv,ncv) , Psize(0), initPsize(0), diffConst1(0), maxSize(0), width(0), height(0), onionList(0) {}
    // N.B. Copy constructor and equals operator do copy the list of onions - much hillarity ensued when they didn't 
    basinObj( const basinObj& b ) : ncv(b.ncv), nonions(b.nonions), det(b.det), R(b.R), refval(b.refval), covar(b.covar), metric(b.metric), Psize(b.Psize), initPsize(b.initPsize), diffConst1(b.diffConst1), maxSize(b.maxSize), width(b.width), height(b.height), onionList(b.onionList) {}
    basinObj& operator=(const basinObj& old){
      if (&old==this) return *this;
      nonions=old.nonions; det=old.det; R=old.R; refval=old.refval; covar=old.covar; metric=old.metric;
      Psize=old.Psize; initPsize=old.initPsize; diffConst1=old.diffConst1; 
      maxSize=old.maxSize; width=old.width; height=old.height; onionList=old.onionList;
      return *this;
    } 
    index size() const { return ncv; } 
    void resize( const index& n ){ ncv=n; refval.resize(ncv); metric.resize(ncv,ncv); covar.resize(ncv,ncv); }
    // Functions used to setup basins
    void setBasin( const index& periodic, const rVector& center, const rMatrix& cov, const real& sizeParam, const real& w, const real& h );
    void setDiffusion( const real& dconst1 ){ diffConst1=dconst1; }
    // Functions for restart
    void restartBasin( const index& periodic, const real& sizeParam, const real& w, const real& h );
    void restartOnion( const onionObj& newOnion, const real& newsize ){ 
      Psize=newsize;
      if( newOnion.hasHeight()==1 ){ onionList.push_back(newOnion); nonions++; return; }
      return;
    }
    // real readOnion( std::istream *onionFile );
    // Functions for metadynamics
    int getNonions() const { return nonions; }
    real calcR( const rVector& colvars, rVector& derivs ) const;
    real calcR( const rVector& period, const rVector& colvars, rVector& derivs ) const;
    real calcPotential( const index& ncv, const rVector& colvars, const rVector& period, rVector& derivatives );
    index inMetadRegion() const {if( R<Psize ){return 1;} return 0;};
    index inExpandRegion() const{if( R>Psize && R<(Psize+width) ){return 1;} return 0;};
    void addOnion( mpiostream *onionfile, const index& bas, const real& time );
    void inflateBasin( mpiostream *onionfile, const real& econst, const index& bas, const real& time );
    // Functions used for diagnostics on cluster analysis
    real getDet() const { return det; }
    real getSigma() const { return Psize/initPsize; }
    // Function used to fit diffusion constants
    real calcD( const rVector& cv1, const rVector cv2 ) const ;
    real calcD( const rVector& period, const rVector& cv1, const rVector cv2 ) const ;
    // Function to return the center of the basin (used by non-linear dimensional reduction algorithms)
    void getCenter( rVector& c ) const { c=refval; }
    // Function used by bespoke colvars only - sets the metric equal to the identity
    void setMetricEqualToIdentity(){ metric=0; for(index i=0;i<ncv;++i){ metric(i,i)=1.0; } }
//#ifdef RECON_DRIVER
//    index inMetadRegion( const rVector& period, const rVector& colvars ) const;
//    index getNHills() const {return nonions;}
//    real calcBias( const real& R ) const;
//    real dist( const rVector& period, const rVector& colvars ) const;
//    real getSize() const { return Psize; }
//#endif
};

std::istream& operator>>(std::istream& istr, basinObj& bas);
std::ostream& operator<<(std::ostream& ostr,const basinObj& bas);

// This calculates the overlap using Mautista's measure
real overlap( const index& ncv, const rVector& period, const basinObj& clus , const basinObj& bas );
// This compares the R values for two basins
index closerBasin( const basinObj& bas1, const basinObj& bas2 );
// This calculates:
//   (1) The distance between bas and clus in the metric of clus
//   (2) The distance between clus and bas in the metric of bas
//   It returns the average of these two distances 
real basinDistance( const index& ncv, const rVector& period, const basinObj& clus , const basinObj& bas );
//#ifdef RECON_DRIVER
//// This calculates the distance between basins
//real dist( const index& ncv, const rVector& period, const basinObj& bas1, const basinObj& bas2 );
//#endif

// This reads in a list of basins
void readBasins( const index& ncv, const index& periodic, const real& sizeParam, const real& width, const real& height, std::vector<basinObj>& basinList );
// This reads in onions
index readOnions( std::vector<basinObj>& basinList );

};  // END OF namespace

#endif

