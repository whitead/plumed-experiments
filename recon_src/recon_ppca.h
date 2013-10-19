#ifndef __rp_clustering
#define __rp_clustering

#include <vector>
#include <valarray>
#include <iostream>
#include <string>
#include <fstream>

#include "recon_types.h"
#include "recon_utils.h"
#include "recon_basins.h"

#define GM_EPS 1E-5
#define MAXTRIES  1000
#define MAXCLUSTERS 10
#define PPCA_SHIFT 0.01
#define PPCA_MIN_COVAR 0.01
#define PPCA_COOLING_RATE 0.95

// We assume this is the periodicity of aperiodic variables when they are mixed with periodic variables.
#define FAKE_PERIOD 100.0

namespace rp {

class ppcaClass;

class ppca_cluster {
     friend std::ostream& operator<<(std::ostream&,const ppca_cluster&);
     friend real distance( const ppca_cluster& , const ppca_cluster& );
     friend class ppcaClass;
private:
     index ncv, dim;
     real weight;
     rVector center;
     rMatrix C;
public:
     ppca_cluster(index ncolvar=0) : ncv(ncolvar), dim(0), weight(0), center(ncv), C(ncv,ncv) {}
     ppca_cluster(const ppca_cluster& c) : ncv(c.ncv), dim(c.dim), weight(c.weight), center(c.center), C(c.C) {}
     ppca_cluster& operator=(const ppca_cluster& old){
       if (&old==this) return *this;
       if( old.ncv!=ncv ) resize( old.ncv );
       dim=old.dim; weight=old.weight; center=old.center; C=old.C; 
       return *this;
     }
     void resize( const index& ncolvar ){ ncv=ncolvar; center.resize(ncv); C.resize(ncv,ncv); }
     void setData( const index& n, const rVector& v, const real& sigma ){ 
        center=v; weight=1./float(n); C=0; for(index i=0;i<ncv;++i){ C(i,i)=sigma; } 
     }
     void getCenter( rVector& c ) const { c=center; } 
     void peturbCenter(){ for(index i=0;i<ncv;++i){ center[i]+=PPCA_SHIFT*( drand48()-0.5 ); } }
     void calcProbs( const index& bas, const rVector& period, const rMatrix& data, rMatrix& probs ) const;
     void calcDAparams( const index& bas, const real& sigma, const index& nclusters, const rVector& period, const rMatrix& data, const rMatrix& probs );
     void calcPPCAparams( const index& bas, const real& sigma, const rVector& period, const rMatrix& data, const rMatrix& probs );
     // This is the number of parameters = number of elements in W + number of elements in center vector + weight of cluster
     index nparams() const { return (dim+1)*ncv+1; }  
     // This returns the fuzzy volume (determinant of the colvar)
     real fuzzyVolume() const { return sqrt( exp(logdet(C)) ); }
     // This returns the weight (used in recon_metad)
     real getWeight() const { return weight; } 
     // This will turn a cluster into a basin
     void extractBasin( const index& periodic_data, const real& sizeParam, const real& width, const real& height, basinObj& basin ) const {
         basin.setBasin( periodic_data, center, C, sizeParam, width, height ); 
     }
};

std::ostream& operator<<(std::ostream& ostr, const ppca_cluster& clust);
real distance( const ppca_cluster& c1, const ppca_cluster& c2 );

class ppcaClass {
private:
    index nCoolSteps;
    index sz, npoints;
    index ncv, lfuzvol;            //, maxclusters;
    rMatrix data, probs;
    rVector period;
    rVector final_weights;
    std::vector<ppca_cluster> clusterList;
    std::vector<ppca_cluster> workingList;
//    std::vector<ppca_cluster> bestList;
public:
    std::ofstream *ofdatafile;
    mpiostream *odatafile;
    ppcaClass() : nCoolSteps(0), sz(0), npoints(0), ncv(0), lfuzvol(0), ofdatafile(NULL), odatafile(NULL) {}
    void setup( const index& ndata, const index& fuzzy, const index& ncolvars, const rVector& per ); // const real& wTol);
    index push_back( const rVector& dat, mpiostream *diagnfile );
    // The Gaussian mixture code
    void gaussMix( mpiostream *diagnfile );
    real getInitSigma() const;
    // void drawNewCenter( const index& ncenters, rVector& newcenter ) const;
    void mStep( const index& rclusters, const index& nclusters, const real& sigma );
    real eStep( const index& nclusters );
    // This extracts everything to the update code
    void getClusters( index& nclus, std::vector<ppca_cluster>& cList ) const { 
       nclus=clusterList.size(); cList.resize(nclus); 
       for(index i=0;i<nclus;++i){ cList[i].resize(ncv); cList[i]=clusterList[i]; }
    }
    // This is used during restart
    index ndatapoints() const { return npoints; }
    // This is used to automatically fit the diffusion constants
    real fitDiffusion( const real& tstep, const rVector& periodicities, const basinObj& bas ) const ;
};

};  //ends namespace rp
    
#endif  // ends #ifndef __rp_clustering
