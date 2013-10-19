#ifdef RECONMETAD 
#include "recon_metad.h"
//#include "recon_ioparser.h"

namespace rp {

// This is used to monitor the basin we are inside during a trajectory
//void reconObj::setupMonitor( const index& ncolvars, const rVector& cv_per ) {
//
//   ncv=ncolvars;     // Number of colective coordinates
//   nbasins=0;        // Used to keep track of the number of basins
//   //tstep=step;       // Store timestep of simulation
//   conserved=0.;     // Conserved quantity
//
//   // Transfer the periodicities of all variables to the local storage space
//   period.resize( ncv ); period=cv_per; index perfix=0;
//   if(period[0]!=0.0){
//       for(index i=1;i<ncv;++i){ if(period[i]==0.0){ perfix=1; break; } } 
//       // for(index i=1;i<ncv;++i){ if(period[i]==0.0) ERROR("YOU CANNOT MIX PERIODIC AND NON PERIODIC VARIABLES"); }
//       // Convert periods to quantities used to calulate the force
//       if(perfix==0){ for(index i=0;i<ncv;++i){ period[i]=2*R_PI/period[i]; } }
//   } else{
//       for(index i=1;i<ncv;++i){ if(period[i]!=0.0){ perfix=1; break; } }
//       // for(index i=1;i<ncv;++i){ if(period[i]!=0.0) ERROR("YOU CANNOT MIX PERIODIC AND NON PERIODIC VARIABLES"); }
//   }
//
//   // Setup fake periodicities 
//   if(perfix==1){ 
//      std::pout<<"IN MIXING PERIODIC AND NON-PERIODIC VARIABLES ASSUMING A PERIOD OF "<<FAKE_PERIOD<<" FOR ALL APERIODIC VARIABLES"<<std::endl;
//      for(index i=0;i<ncv;++i){ if(period[i]==0.0){ period[i]=FAKE_PERIOD; } period[i]=2*R_PI/period[i]; }
//   }
//
//   // Read in the basins data
//   index periodic_data=0; if(period[0]!=0.0){periodic_data=1;} 
//   readBasins( ncv, periodic_data, sizeParam, width, height, basinList );
//   nbasins=basinList.size();
//
//   // Read in the onions data
//   conserved=readOnions( basinList ); 
//}
//
//int reconObj::doMonitor( const rVector& colvars ){
//
//   index flag=0, curbas; real minvol=0, vol; rVector derivatives( ncv );
//   // std::cout<<"I AM INSIDE THE FOLLOWING BASINS ";
//   for(index i=0;i<nbasins;++i){
//      // We must compute R first
//      basinList[i].calcPotential( ncv, colvars, period, derivatives );
//
//      if( basinList[i].inMetadRegion()==1 && flag==0 ){ flag=1; curbas=i; minvol=ncv*log( basinList[i].getSigma() )+basinList[i].getDet(); }   //std::cout<<curbas<<" "; }
//      else if( basinList[i].inMetadRegion()==1 ){
//          //std::cout<<i<<" ";
//          vol=ncv*log( basinList[i].getSigma() )+basinList[i].getDet();
//          if( vol<minvol ){ minvol=vol; curbas=i; }
//      }
//   }
//   //std::cout<<std::endl;
//   if( flag==1 ){ return curbas; }
//   return -1;
//}

void reconObj::setup( const index& ncolvars, const rVector& cv_per, const real& step, struct recon_data_s& reconinpt, FILE* fplog ){

   srand48(12345);   // Initialize the random number generator
   ncv=ncolvars;     // Number of colective coordinates
   nbasins=0;        // Used to keep track of the number of basins
   tstep=step;       // Store timestep of simulation

   // Transfer the periodicities of all variables to the local storage space
   period.resize( ncv ); period=cv_per; index perfix=0;
   if(period[0]!=0.0){
       for(index i=1;i<ncv;++i){ if(period[i]==0.0){ perfix=1; break; } }
       // for(index i=1;i<ncv;++i){ if(period[i]==0.0) ERROR("YOU CANNOT MIX PERIODIC AND NON PERIODIC VARIABLES"); }
   } else{
       for(index i=1;i<ncv;++i){ if(period[i]!=0.0){ perfix=1; break; } } 
       // for(index i=1;i<ncv;++i){ if(period[i]!=0.0) ERROR("YOU CANNOT MIX PERIODIC AND NON PERIODIC VARIABLES"); }
   }

   // Setup fake periodicities
   if(perfix==1){
      WARNING("IN MIXING PERIODIC AND NON-PERIODIC VARIABLES ASSUMING A PERIOD OF "<<FAKE_PERIOD<<" FOR ALL APERIODIC VARIABLES"); 
      for(index i=0;i<ncv;++i){ if(period[i]==0.0){ period[i]=FAKE_PERIOD; } }
   }

   index restartFlag; restartFlag=reconinpt.restart;    // Flag for restart
   // Onion parameters
   height=reconinpt.height;      // Height of onions
   stride=reconinpt.stride;      // Frequency of generation of onions
   width=reconinpt.width;        // Width of onions
   //temperature=reconinpt.temp;   // System temperature

   // Check everything has been set
   if( height<0.0 ) ERROR("NO ONION HEIGHT SET");
   if( stride<0 ) ERROR("NO ONION STRIDE SET");
   if( width<0.0 ) ERROR("NO ONION WIDTH SET");
//   if( temperature<0.0 ) ERROR("SYSTEM TEMPERERATURE IS NOT SET");

   //std::pout<<"INPUT FOR RECONNAISSANCE METADYNAMICS"<<std::endl;
   //std::pout<<"METADYNAMICS PARAMETERS: HEIGHT "<<height<<" WIDTH "<<width<<" STRIDE "<<stride<<std::endl;
   //<<" SYSTEM TEMPERATURE "<<temperature<<std::endl;

   // Basin parameters
   basTol=reconinpt.basTol;               // Threshold for accepting basins
   sizeParam=reconinpt.iSize;             // Parameter for setting initial size of basins
   expandParam=reconinpt.diffconst;       // This controls the rate of expansion of the basins         
//   expandParam=( reconinpt.diffconst*stride*tstep) / (2*width);   // This controls the rate of expansion of the basins

   //if(expandParam>0){
   //   std::pout<<"BASIN UPDATE PARAMETERS: TOLERANCE "<<basTol<<" INITIAL SIZE "<<sizeParam<<" EXPAND PARAM "<<reconinpt.diffconst<<std::endl;
   //} else {
   //   std::pout<<"BASIN UPDATE PARAMETERS: TOLERANCE "<<basTol<<" INITIAL SIZE "<<sizeParam<<" FITTING DIFFUSION CONSTANT"<<std::endl;
   //}

   // Check everythign has been set
   if( basTol<0.0 ) ERROR("TOLERANCE FOR ACCEPTING BASINS NOT SET");
   if( sizeParam<0.0 ) ERROR("INITIAL SIZE OF BASINS NOT SET");
   if( expandParam<0.0 ) ERROR("EXPANSION PARAMETER NOT SET");

   // Now set up clustering objects
   nclusterObj=reconinpt.nscales; myclusteringObjs.resize( nclusterObj ); myclusterFreqs.resize( nclusterObj );
   //std::pout<<"DOING CLUSTER ANALYSIS AT "<<nclusterObj<<" TIME SCALES"<<std::endl;
   for(index i=0;i<nclusterObj;++i){
      if( reconinpt.runFreq[i]<0 ) ERROR("FREQUENCY TO RUN CLUSTER ANALYSIS NOT SET "<<i<<"th CLUSTER OBJECT");
      if( reconinpt.storeFreq[i]<0 ) ERROR("FREQUENCY TO STORE CLUSTERS NOT SET "<<i<<"th CLUSTER OBJECT");
      if( reconinpt.runFreq[i]%reconinpt.storeFreq[i]!=0 ) ERROR("CLUSTER ANALYSIS FREQUENCY MUST BE A MULTIPLE OF THE CLUSTER STORAGE FREQUENCY "<<i<<"th CLUSTER OBJECT");
      // Setup the clustering object
      reconinpt.runFreq[i]/=reconinpt.storeFreq[i]; myclusterFreqs[i] = reconinpt.storeFreq[i];
   //   std::pout<<"CLUSTER ANALYSIS SCALE "<<i+1<<" PARAMETERS: STORE FREQUENCY "<<myclusterFreqs[i]<<" ";
      // This is the old pca + gm
//      myclusteringObjs[i].setup( reconinpt.runFreq[i], ncv, period, reconinpt.clusterTol[i], reconinpt.ntries[i], reconinpt.pcaTol[i], basTol );  
      // This is the new mixture of ppca models
      myclusteringObjs[i].setup( reconinpt.runFreq[i], reconinpt.fuzzy, ncv, period );
   }

   // Output parameters for reconnaissance metadynamics
   fprintf(fplog, "|--- RECONNAISSANCE METADYNAMICS PARAMETERS --- \n \n");
   fprintf(fplog, "|-DOING CLUSTER ANALYSIS AT %d TIME SCALES\n",nclusterObj);
   for(index i=0;i<nclusterObj;++i){
     fprintf(fplog, "|-CLUSTERING TIMESCALE %d : STORE FREQUENCY %d RUN FREQUENCY %d\n",i+1,reconinpt.storeFreq[i], reconinpt.runFreq[i]);
   }
   fprintf(fplog,"\n");
   fprintf(fplog,"|-ONIONS : WIDTH %f, HEIGHT %f, WRITING STRIDE %d \n", width, height, stride);
   fprintf(fplog,"|-BASINS : TOLERANCE %f, INITIAL SIZE %f, EXPAND PARAM %f \n", basTol, sizeParam, expandParam);

   // Now convert periods to the quantity used when calculating the force 
   if(period[0]!=0.0){
      for(index i=0;i<ncv;++i) period[i]=2*R_PI/period[i];
   }

   if(restartFlag==1){
      // Call restart
      index nonions; nonions=restart();
      // Open the CLUSTER_DATA files in append mode
      for(index i=0;i<nclusterObj;++i){
         myclusteringObjs[i].ofdatafile=new std::ofstream(
                     (std::string("CLUSTER_DATA.")+int2str(i)).c_str(),
                     std::ios_base::app);
         myclusteringObjs[i].odatafile=new mpiostream(*(myclusteringObjs[i].ofdatafile));
      }
      fprintf(fplog,"\n");
      fprintf(fplog,"|-DURING RESTART READ IN %d BASINS, %d ONIONS \n",basinList.size(),nonions);
      for(index i=0;i<nclusterObj;++i){
         fprintf(fplog,"|-RESTARTING %dth CLUSTERING TIMESCALE WITH %d POINTS READ FROM CLUSTER.DATA.%d \n",i+1, myclusteringObjs[i].ndatapoints(), i);
      }
   }
   else{
      // Open the various output files and write over any old ones
      basinfile=new mpiostream(*(new std::ofstream("BASINS")));
      onionfile=new mpiostream(*(new std::ofstream("ONIONS")));
      diagnfile=new mpiostream(*(new std::ofstream("PPCA_DIAGNOSTICS"))); 
      for(index i=0;i<nclusterObj;++i){
         myclusteringObjs[i].ofdatafile=new std::ofstream( (std::string("CLUSTER_DATA.")+int2str(i)).c_str() );
         myclusteringObjs[i].odatafile=new mpiostream(*(myclusteringObjs[i].ofdatafile)); 
      }
   }
   fprintf(fplog,"\n");

   // Tempory measure to write out forces
   // forcefile=new mpiostream(*(new std::ofstream("FORCES")));
}

real reconObj::doReconMetad( const int timestep, const rVector& colvars, rVector& derivatives ){

   // This is everything for cluster analysis
   for(index i=0;i<nclusterObj;++i){  
      if(timestep%myclusterFreqs[i]==0){
         if(myclusteringObjs[i].push_back( colvars, diagnfile )==0){ 
            real time; time=timestep*tstep;
            updateBasinList(i,time); 
            myclusteringObjs[i].ofdatafile->close();
            delete myclusteringObjs[i].odatafile;
            delete myclusteringObjs[i].ofdatafile;  
            myclusteringObjs[i].ofdatafile=new std::ofstream(
                     (std::string("CLUSTER_DATA.")+int2str(i)).c_str() );  
            // myclusteringObjs[i].ofdatafile=new std::ofstream("CLUSTER_DATA");
            myclusteringObjs[i].odatafile=new mpiostream(*myclusteringObjs[i].ofdatafile);  
         }
      }
   }

   // This calculates the potential - N.B. only those hills with onions need be considered
   real ene=0.0; derivatives=0.0;
   for(index i=0;i<nbasins;++i){
       if(basinList[i].getNonions()>0) ene+=basinList[i].calcPotential( ncv, colvars, period, derivatives );
   }

   // (*forcefile)<<timestep*tstep<<" "<<derivatives<<std::endl;

   // This adds onions
   if( timestep%stride==0 ){
      real time,tmp; time=timestep*tstep; rVector fake_der(ncv); index nRegions=0, nExpand=0, bas;
      for(index i=0;i<nbasins;++i){
         // This ensures that the R for hills with no onions gets calculated
         if(basinList[i].getNonions()==0) tmp=basinList[i].calcPotential( ncv, colvars, period, fake_der );
         // This is a counter over how many basins we are inside (we add onions to the closest of these)
         if( basinList[i].inMetadRegion()==1 ){ 
            nRegions++;
            if( nRegions==1 ){ bas=i; }else if(closerBasin( basinList[i], basinList[bas] )==1){ bas=i; } 
         }
         // This is a check if there are any basins we can expand (only do this if we are in no basins)
         else if( nRegions==0 && basinList[i].inExpandRegion()==1 ){ nExpand++; bas=i; }
      }
      // I am phasing out the conserved quantity as we never used it
      // if(nRegions>=1){ conserved+=basinList[bas].addOnion( onionfile, bas, time ); }

      // This does the inflation of basins and the addition of onions
      if(nRegions>=1){ basinList[bas].addOnion( onionfile, bas, time ); }
      // We expand one basin at a time
      else if(nExpand==1){ 
         // Calculate the number we must multiply diffusion constants by to ge the expansion parameter 
         real expandConst; expandConst = ( stride * tstep ) / ( 2*width ); 
         basinList[bas].inflateBasin( onionfile, expandConst,  bas, time ); 
      }
   } 

   return ene;
}

void reconObj::updateBasinList( const index& objno, const real& time ){
   real overl, maxoverlap, thisoverlap; basinObj newbasin( ncv );

   // Get the real periodicities rather than the factors used for calculating the force
   rVector periodicities( ncv ); periodicities=0.; index periodic_data=0;
   if(period[0]!=0.0){
     periodic_data=1;
     for(index i=0;i<ncv;++i) periodicities[i]=2*R_PI/period[i];
   }

   // std::pout<<"RESULTS FOR PPCA ANALYSIS AT TIME "<<time<<std::endl;

   // Get the basins from the clustering algoirthm
   index ntmpb; std::vector<ppca_cluster> clusterList; 
   myclusteringObjs[objno].getClusters( ntmpb, clusterList );

   // Now work through basins and establish which basins to add
   (*diagnfile)<<"FATES OF FITTED BASINS"<<std::endl;
   for(index i=0;i<ntmpb;++i){
      overl=clusterList[i].getWeight();
      if( overl<basTol ){ (*diagnfile)<<" OVERLAP "<<overl<<" REJECTED"<<std::endl; continue; }
      // Extract a basin from the list of clusters
      clusterList[i].extractBasin( periodic_data, sizeParam, width, height, newbasin );
      // Calculate the overlap with the old basins
      maxoverlap=0.0;
      for(index j=0;j<basinList.size();++j){
         thisoverlap=overlap( ncv, periodicities, newbasin, basinList[j] );
         if( thisoverlap>maxoverlap ) maxoverlap=thisoverlap;
      }
      overl*=1.0-maxoverlap;
      (*diagnfile)<<" OVERLAP  "<<overl;
      // Add to basin List if appropriate
      if( overl>basTol ){
          (*diagnfile)<<" ACCEPTED AS BASIN "<<basinList.size()<<" "; 
          // if( expandParam>0 ){ newbasin.setDiffusion( expandParam ); std::pout<<"\n"; }
          // else { 
          // Calculate diffusion constant from clustering data
          real diffConst; diffConst=expandParam*myclusteringObjs[objno].fitDiffusion( (myclusterFreqs[objno]*tstep) , period, newbasin );
          //real initSize; initSize=sqrt( ncv-1.0 ) + sizeParam; 
          //real dconst; dconst = expandParam - (diffConst*stride*tstep) / ( 2*width*initSize );
          (*diagnfile)<<"DIFFUSION CONSTANT IN BASIN EQUALS "<<diffConst<<std::endl; 
          // Set diffusion constant in basin
          newbasin.setDiffusion( diffConst );                         // , dconst ); 
          // }
          basinList.push_back( newbasin );   // Add to list of basins
          (*basinfile)<<basinList.size()<<" "<<time<<" "<<newbasin<<std::endl;  // Write to basins file  
      }
      else{(*diagnfile)<<" REJECTED"<<std::endl;}
   }

   (*diagnfile)<<"====END OF PPCA DIAGNOSTICS FOR CLUSTERING AT TIME "<<time<<"====\n \n";
   (*diagnfile).flush();

   nbasins=basinList.size();      // Update the number of basins
}

index reconObj::restart(){
   // Read in the basins
   index periodic_data=0; if(period[0]!=0.0){periodic_data=1;}
   readBasins( ncv, periodic_data, sizeParam, width, height, basinList );
   nbasins=basinList.size();

   // Read in the onions and store the conserved quantity
   index nonions; nonions=readOnions( basinList );

   // Reopen the output files
   basinfile=new mpiostream(*(new std::ofstream("BASINS",std::ios_base::app)));
   onionfile=new mpiostream(*(new std::ofstream("ONIONS",std::ios_base::app)));
   diagnfile=new mpiostream(*(new std::ofstream("PPCA_DIAGNOSTICS",std::ios_base::app)));

   // Read in the cluster data files 
   rVector dat(ncv); index tstep; std::string line; real extranum;
   for(index i=0;i<nclusterObj;++i){
      std::ifstream* idatafile=new std::ifstream( (std::string("CLUSTER_DATA.")+int2str(i)).c_str() );
      if (!idatafile->good()) ERROR("Something bad happened while opening "<<i<<"th CLUSTER_DATA file for restart");
      while( getline( (*idatafile), line) ) {
         // Put the line we have read in onto a string stream and read in
         std::stringstream input (line); // input>>tstep>>dat;
         if( ( input >> tstep >> dat ).fail() ) ERROR("WRONG FORMAT IN CLUSTER DATA FILE");
         // This should test that what we read is sensible
         if( !(input>>extranum).fail() ) WARNING("STRANGE NUMBER OF COLUMNS IN CLUSTER DATA FILE") 
         // Push back coordinates to the cluster object 
         if( myclusteringObjs[i].push_back( dat , diagnfile )==0 ){ updateBasinList(i,-1000.); }
      }
      // while( idatafile->good() ) {
      //    (*idatafile)>>tstep>>dat;
      //    if ( !idatafile->good() ) break;
      //    if( myclusteringObjs[i].push_back(dat)==0 ){ updateBasinList(i,-1000.); } 
      // }
//      std::pout<<"Read in "<<myclusteringObjs[i].ndatapoints()<<" data points from "<<i<<"th CLUSTER_DATA file"<<std::endl;
      idatafile->close(); delete idatafile; 
   }
   return nonions;
}

}  // END OF NAMESPACE
#endif
