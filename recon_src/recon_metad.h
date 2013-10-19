#ifndef __rp_recon
#define __rp_recon

#include <stdio.h>
#include <vector>
#include <valarray>
#include <iostream>
#include <string>
#include <fstream>

#include "recon_types.h"
#include "recon_utils.h"
// #include "recon_metric.h"
// #include "recon_cluster_gmm.h"
// #include "recon_cluster_hgmm.h"
#include "recon_ppca.h"
#include "recon_basins.h"

extern "C" {
#include "recon_cbind.h"                       
}

namespace rp {

// #define cluster gmm_cluster
// #define clusteringClass gmmClass

class reconObj{
private:
    // Output files
    mpiostream* basinfile;
    mpiostream* onionfile; 
    mpiostream* diagnfile;
    // Forces output
    // mpiostream* forcefile;
    // Conserved quantity
    // real conserved;
    // Parameters
    rVector period;
    index ncv, nbasins;
    real basTol, tstep; //, boltz, Qpadding;
    real sizeParam, expandParam;
    index stride;
    real height, width; //  temperature; 
    std::vector<basinObj> basinList;
    index nclusterObj;
    std::vector<index> myclusterFreqs;
//    std::vector<clusteringClass> myclusteringObjs;
    std::vector<ppcaClass> myclusteringObjs;
public:
    reconObj() : ncv(0), nbasins(0), basTol(0), tstep(0), sizeParam(0), expandParam(0), stride(0), height(0), width(0), nclusterObj(0) {}
    void setup( const index& ncolvars, const rVector& cv_per, const real& step, struct recon_data_s& reconinpt, FILE* fplog );
    real doReconMetad( const int timestep, const rVector& colvars, rVector& derivatives );
//    real getConservedQuantity(){return conserved;}
    void updateBasinList(const index& objno, const real& time);
    index restart();
    // These are used for testing derivatives only
    void setup_test_derivatives( const rVector& colvars );
    void test_derivatives_calc_cv( const rVector& colvars, real& pos );
    void test_derivatives_calc_derivs( const rVector& colvars, rVector& derivs );
    // This is for monitoring what basins we enter during unbiassed MD
    //void setupMonitor( const index& ncolvars, const rVector& cv_per );
    //int doMonitor( const rVector& colvars );

    // These are used for testing the onion derivatives  (I can't remember what these did so I got rid of them
//    void setup_test_onions( const rVector& colvars );
//    void test_onions( rVector& colvars ); 
};

};  // End of namespace

#endif  // ends #ifndef __rp_recon
