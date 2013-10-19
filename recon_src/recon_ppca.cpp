#ifdef RECONMETAD 
#include "recon_ppca.h"

namespace rp {

std::ostream& operator<<(std::ostream& ostr,const ppca_cluster& clust){
   ostr<<" dim "<<clust.dim<<" center "<<clust.center<<" covar "<<clust.C<<" weight "<<clust.weight; // For time being just write position
   return ostr;
}

real distance( const ppca_cluster& c1, const ppca_cluster& c2 ){
#ifdef DEBUG
  if( c1.ncv!=c2.ncv ) ERROR("clusters do not have same dimensionality");
#endif

  real sum=0., tmp;
  for(index i=0;i<c1.ncv;++i){ tmp=c1.center[i]-c2.center[i]; sum+=tmp*tmp; }
  return sqrt(sum);

}

void ppcaClass::setup( const index& ndata, const index& fuzzy, const index& ncolvars, const rVector& per ){      //, const real& wTol ){
   
   npoints=0;                  // Counter on the number of points in the data matrix
   sz=ndata;                   // When we have this number of points do cluster analysis
   ncv=ncolvars;               // Total number of collective coordinates we are employing
   lfuzvol=fuzzy;              // Are we using fuzzy volume to select the number of clusters
   period.resize( ncv ); period=per; // Copy the periodicities of the variables to the class variable here

   // Now check for fake periods and remove them
   for(index i=0;i<ncv;++i){ 
      if(period[i]==FAKE_PERIOD){ period[i]=0.0; }
   }

   data.resize( ndata, ncv );             // Resize data array 
   
   // Now get the stuff from input 
   // nTries=tries;                          // Number of tries at each size
   nCoolSteps=int( log(PPCA_MIN_COVAR) / log(PPCA_COOLING_RATE) ) + 1;
   // minWeight=wTol;                        // Minimum weight any cluster should have
   // pdim=dim;                              // Number of principal compnents used in fitting the covariance 
   // maxclusters=int(1.0/wTol)+1;      // Maximum number of clusters possible with this minimum weight
   
   // std::pout<<"NSNAPS "<<sz<<" COOLING FOR "<<nCoolSteps<<" STEPS"<<std::endl;

   probs.resize( MAXCLUSTERS,ndata );     // Setup a probabilities array
   workingList.resize( MAXCLUSTERS );     // This stores clusters during cluster algorithm
//   bestList.resize( maxclusters );        // This stores the best clusters 
   // Set the correct number of colvars for all the clusters
   for(index n=0;n<MAXCLUSTERS;++n){ workingList[n].resize(ncv); }
}

index ppcaClass::push_back( const rVector& dat , mpiostream *diagnfile ){
   // Store the vector of collective coordinates
   data.setRow(npoints,dat); npoints++;
   // Write out points in cluster data file (so we can restart properly if we crash)
   if (odatafile!=NULL) { (*odatafile)<<npoints<<" "<<dat<<std::endl; }
   // If we have a full datamatrix do a cluster analysis
   if(npoints==sz){ gaussMix( diagnfile ); npoints=0; }
   // Otherwise print the data so that restart will work
   // else if (odatafile!=NULL) { (*odatafile)<<npoints<<" "<<dat<<std::endl; }
   // Return the number of points
   return npoints;
}

void ppcaClass::gaussMix( mpiostream *diagnfile ){
   real oldloglike, loglike, sigma0, sigma, bik, oldbik, fvol;
   rVector newcenter( ncv ); index flag;

   // Calculate maximum eigenvector of the covariance
   // This is used to set the intial temperature
   sigma0=getInitSigma();

   // Set the probability of all basins equal to one
   for(index i=0;i<sz;++i){ probs[0,i]=1.0;}
   // Resize the cluster list 
   clusterList.resize(1); clusterList[0].resize(ncv);

   // Do a PPCA fit of the data with one basin only N.B. first 1 equals do PPCA
   mStep(1, 1,  sigma0*PPCA_MIN_COVAR ); clusterList[0]=workingList[0];
   // Compute the bayesian information criterion
   oldbik=clusterList[0].nparams()*log( sz ) - 2*eStep( 1 );
   fvol=clusterList[0].fuzzyVolume();

   if(lfuzvol==1){ oldbik=fvol; }

   //std::pout<<"CLUSTER ANALYSIS DIAGNOSTICS"<<std::endl;
   (*diagnfile)<<"COMPLETED FIT WITH 1 CLUSTER. BIK EQUALS "<<oldbik<<" FUZZY VOLUME "<<fvol<<std::endl; 

   // Get position of center of cluster
   clusterList[0].getCenter( newcenter );

   // Loop over number of clusters
   for(index nclusters=2;nclusters<=MAXCLUSTERS;++nclusters){
       sigma=sigma0; flag=0; // Set initial temperature
       // Put all clusters at center at start
       for(index m=0;m<nclusters;++m) workingList[m].setData( nclusters, newcenter, sigma ); 

       // Now begin the simulated annealing
       for(index k=0;k<nCoolSteps;++k){
          // Peturb all the clusters if we are doing DA
          if(flag==0){ for(index m=0;m<nclusters;++m) workingList[m].peturbCenter(); }
          // Do a first eStep
          oldloglike=eStep( nclusters );  
          // This actually does expectation maximization 
          for(index kk=0;kk<MAXTRIES;++kk){
             mStep( flag, nclusters , sigma );
             loglike=eStep( nclusters );
             if( fabs( (loglike-oldloglike)/loglike ) < GM_EPS ){ break; }
             oldloglike=loglike;
          }
          sigma*=PPCA_COOLING_RATE; 

          // Check for end of DA 
          if(flag==0){
             flag=1;
             for(index m=1;m<nclusters;++m){
                for(index n=0;n<m;++n){
                    if(distance( workingList[m],workingList[n] )<ncv*PPCA_SHIFT) flag=0;
                }
             }
             //if(flag==1) std::pout<<"SWITCHING TO PPCA AT STEP "<<k<<" SIGMA EQUALS "<<sigma<<std::endl;
          }
       }
       // Compute the Bayesian Information Criterion
       bik=0.;for(index m=0;m<nclusters;++m){ bik+=workingList[m].nparams(); } bik*=log(sz); bik-=2*loglike;
       fvol=0;for(index m=0;m<nclusters;++m){ fvol+=workingList[m].fuzzyVolume(); }
       if( flag==0 ){ (*diagnfile)<<"FAILED TO SEPARATE CLUSTERS DURING DA ANNEALING : BIK FOR ISOTROPIC MODEL "<<bik<<" FUZZY VOLUME "<<fvol<<std::endl; }
       else{ (*diagnfile)<<"COMPLETED FIT WITH "<<nclusters<<" CLUSTERS. BIK EQUALS "<<bik<<" FUZZY VOLUME "<<fvol<<std::endl; }

       if(lfuzvol==1){ bik=fvol; }      

       // We want to take the model with the maximum bik (in which clusters have split up)
       if( bik<oldbik && flag!=0 ){
           oldbik=bik; clusterList.resize(nclusters); 
           for(index m=0;m<nclusters;++m){ clusterList[m].resize(ncv); clusterList[m]=workingList[m]; }
       }
   }
   if(clusterList.size()==MAXCLUSTERS) WARNING("BIC MAXIMUM WAS AT MAXCLUSTERS, CONSIDER INCREASING MAXCLUSTERS IN CODE");
   (*diagnfile)<<"BEST FIT FOR DATA USES "<<clusterList.size()<<" GAUSSSIANS"<<std::endl;
//   for(index m=0;m<clusterList.size();++m){ std::pout<<"CLUSTER "<<m<<" "<<clusterList[m]<<std::endl; }

}

real ppcaClass::getInitSigma( ) const {

   rVector center(ncv); real sum, sumdiag=0; rMatrix covar(ncv,ncv);
   for(index i=0;i<ncv;++i){
      // This does the average for non periodic variables
      if(period[i]==0){
         sum=0.; 
         for(index k=0;k<sz;++k) sum+=data(k,i);
         center[i]=sum/float(sz);
      }    
      // This does the averae for periodic variables
      else {  
         real phi,sinsum,cossum,pfactor; sinsum=0.0; cossum=0.0; pfactor=2*R_PI / period[i];
         for(index k=0;k<sz;++k){ phi=data(k,i)*pfactor; sinsum+=sin(phi); cossum+=cos(phi); }
         center[i]=atan2( sinsum/float(sz), cossum/float(sz) ) / pfactor;
      }          
      // This does the covariance
      for(index j=0;j<=i;++j){
          sum=0.0; for(index k=0;k<sz;++k) sum+=pbc((data(k,i)-center[i]),period[i])*pbc((data(k,j)-center[j]),period[j]);
          covar(j,i)=covar(i,j)=sum/float(sz);
      }
   }

   // Now diagonalize the covariance to get eigenvalues and eigenvectors
   rVector eigenvals( ncv); rMatrix eigenvecs( ncv,ncv );
   diagMat( covar, eigenvals, eigenvecs );

   // Get the maximum eigenvalue
   return eigenvals[ncv-1];
}

void ppcaClass::mStep( const index& stage, const index& nclusters, const real& sigma ){
   // if the clusters are not separated then we use DA
   if( stage==0 ){
      for(index i=0;i<nclusters;++i){ workingList[i].calcDAparams( i, sigma, nclusters, period, data, probs); }
   }
   // once all clusters are separated calculate using the PPCA tool
   else {
      for (index i=0;i<nclusters;++i){ workingList[i].calcPPCAparams( i, sigma, period, data, probs); }
   }
}

void ppca_cluster::calcDAparams( const index& bas, const real& sigma, const index& nclusters, const rVector& period, const rMatrix& data, const rMatrix& probs ){
   real sum; index np; np=data.rows();

   // First sum up the total probability of this cluster 
   real wgt=0;
   for(index j=0;j<np;++j) wgt+=probs(bas,j);
//   // This computes the fraction of the distribution that is in this cluster 
//   weight=wgt/float(np);

   // In DA weights of all clusters should be the same
   weight = 1.0 / float(nclusters); 

   // Compute the position of the center
   for(index i=0;i<ncv;++i){
      // This does the average for non periodic variables
      if(period[i]==0){
        sum=0.; for(index k=0;k<np;++k) sum+=probs(bas,k)*data(k,i);
        center[i]=sum/wgt;
      }
      // This does the average for periodic variables
      else {
        real phi,sinsum,cossum,pfactor; sinsum=0.0; cossum=0.0; pfactor=2*R_PI / period[i];
        for(index k=0;k<np;++k){ phi=data(k,i)*pfactor; sinsum+=probs(bas,k)*sin(phi); cossum+=probs(bas,k)*cos(phi); }
        center[i]=atan2(sinsum/wgt, cossum/wgt) / pfactor;
      }
   }

   // And set the covariance so that the probability distribution is spherical 
   C=0; for(index i=0;i<ncv;++i){ C(i,i)=sigma; } dim=0;
}

void ppca_cluster::calcPPCAparams( const index& bas, const real& sigma, const rVector& period, const rMatrix& data, const rMatrix& probs ){
   real sum; index np; np=data.rows();  // Number of data points
   rMatrix covar( ncv, ncv ), eigenvecs( ncv, ncv); rVector eigenvals( ncv ); 

   // First sum up the total probability of this cluster 
   real wgt=0; 
   for(index j=0;j<np;++j) wgt+=probs(bas,j);
   // This computes the fraction of the distribution that is in this cluster 
   weight=wgt/float(np);

   for(index i=0;i<ncv;++i){
      // This does the average for non periodic variables
      if(period[i]==0){
         sum=0.;
         for(index k=0;k<np;++k) sum+=probs(bas,k)*data(k,i); 
         center[i]=sum/wgt;
      }
      // This does the average for periodic variables
      else {
         real phi,sinsum,cossum,pfactor; sinsum=0.0; cossum=0.0; pfactor=2*R_PI / period[i];
         for(index k=0;k<np;++k){ phi=data(k,i)*pfactor; sinsum+=probs(bas,k)*sin(phi); cossum+=probs(bas,k)*cos(phi); }
         center[i]=atan2(sinsum/wgt, cossum/wgt) / pfactor;
      }
      // This does the covariance
      for(index j=0;j<=i;++j){
         sum=0.0; for(index k=0;k<np;++k) sum+=pbc((data(k,i)-center[i]),period[i])*pbc((data(k,j)-center[j]),period[j])*probs(bas,k);
         covar(j,i)=covar(i,j)=sum/wgt;
      }
   }

   // Now diagonalize the covariance to get eigenvalues and eigenvectors
   diagMat( covar, eigenvals, eigenvecs );

   // Compute the dimensionality (number of eigenvalues greater than sigma)
   index n; dim=0;
   for(index j=0;j<ncv;++j){ 
      n=ncv-1-j; 
      if( eigenvals[n]>sigma ){ dim++; }
      else{ break; }
   }

   // Set everything to the correct size
   rMatrix U(ncv,dim), W(ncv,dim); rVector Lambda(dim); 
   // Now calculate W
   C=0.;
   for(index j=0;j<dim;++j){
      n=ncv-1-j; Lambda[j]=sqrt(eigenvals[n]-sigma);
      for(index k=0;k<ncv;++k) U(k,j)=eigenvecs(n,k);
   }

   // Complete calculation of W and calculate C
   if(dim>0){ matDMatMult( U, Lambda, W ); matTMatMult( W, W, C ); }
   // Add sigma if we are not in full dimensionality space
   for(index j=0;j<ncv;++j) C(j,j)+=sigma;
}

real ppcaClass::eStep( const index& nclusters ){

   for(index i=0;i<nclusters;++i){ workingList[i].calcProbs( i, period, data, probs ); } 

   real loglike=0. ,sum, tmp, max;
   for(index i=0;i<sz;++i){
     // Find the maximum "probability" in the data
     max=-99e99; for(index j=0;j<nclusters;++j){ if ( probs(j,i) > max ) max=probs(j,i);}
     // Sum the points using the log-sum-exp formula ( log[ \sum exp(z) ] = z_max + log[ \sum exp(z - z_max) ] 
     sum=0.0; for(index j=0;j<nclusters;++j) sum+=exp( probs(j,i) - max );
     tmp = max + log(sum);
     // Normalize by deviding by total probability (the exponent of the log sum) 
     for(index j=0;j<nclusters;++j){
        probs(j,i)-=tmp; probs(j,i)=exp( probs(j,i) );
     }
     loglike+=tmp;  // Add the log of the probability for being at this point to the log likelihood
   }
     
   return loglike;
}

void ppca_cluster::calcProbs( const index& bas, const rVector& period, const rMatrix& data, rMatrix& probs ) const {
   // Quick way to get the logarithm of the determinant of covar (used in normalization)
   rMatrix chol; cholesky( C, chol );
   real lndet=0; for(index i=0;i<ncv;++i){ lndet+=log( chol(i,i) ); } lndet*=2.0;

   rVector u( ncv ), v( ncv ); real loglike=0.;
   for(index i=0;i<data.rows();++i){
      // Compute the distance between this point and the center
      for(index j=0;j<ncv;++j) u[j] = pbc(data(i,j) - center[j], period[j]);
      // Now use cholesky trick to compute u[m] times the inverse covariance
      chol_elsolve( chol, u, v );
      // And compute the unormalizedcluster probability of this point in the distribution
      probs(bas,i)= -0.5*( dotProduct(v,v) + lndet + ncv*log(2.0*R_PI) ) + log( weight );
   } 
}

// This fits our diffusion constants
real ppcaClass::fitDiffusion( const real& tstep, const rVector& periodicities, const basinObj& bas ) const {
   real D, tmp; rVector cv1(ncv), cv2(ncv);

   D=0.0;
   for(index i=1;i<data.rows();++i){
      data.getRow( i-1, cv1 ); data.getRow( i, cv2 );
      if(periodicities[0]!=0){ tmp=bas.calcD( periodicities, cv1, cv2 ); }
      else{ tmp=bas.calcD( cv1, cv2 ); }
      D+=tmp;
   }
   return D / (tstep*(data.rows()-1)*ncv); 
} 

}   // End of namespace rp

#endif
