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
#include "metadyn.h"
#include <assert.h>

void PREFIX read_restraint(struct mtd_data_s *mtd_data)
{
  double uno, due, tre, quattro;
  int i, j, icv, count, iw, nw, ix, iline, tmpc;
  FILE *file;
  char metafile[120], **word, tmpmeta[120];

// object containing parsed input
  t_plumed_input input;

// open PluMeD parameters file  
  file = fopen(mtd_data->metaFilename, "r");
  fflush(mtd_data->fplog);
  if(file==NULL) {
#if defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
    if( (mtd_data->repl==-1)&& (mtd_data->mcr->ms)==0 ){
      char buf[1024];
      sprintf(buf, "MISSING PLUMED INPUT FILE %s",mtd_data->metaFilename);
      plumed_error(buf);
    } else if(mtd_data->mcr->ms!=0 &&mtd_data->repl==-1) {
      strcpy(tmpmeta,mtd_data->metaFilename);
      tmpmeta[strlen(mtd_data->metaFilename) - 4] = '\0';
      sprintf(tmpmeta+strlen(tmpmeta),"%d",mtd_data->mcr->ms->sim);
      sprintf(metafile, "%s.dat", tmpmeta);
      file = fopen(metafile, "r");
      if(file==NULL) {
        char buf[1024];
        sprintf(buf, "MISSING PLUMED INPUT FILE %s",mtd_data->metaFilename);
        plumed_error(buf);
      }
    } else if(mtd_data->repl>-1) {
      strcpy(tmpmeta,mtd_data->metaFilename);
      tmpmeta[strlen(mtd_data->metaFilename) - 4] = '\0';
      sprintf(tmpmeta+strlen(tmpmeta),"%d",mtd_data->repl);
      sprintf(metafile, "%s.dat", tmpmeta);
      file = fopen(metafile, "r");
      if(file==NULL) {
        char buf[1024];
        sprintf(buf, "MISSING PLUMED INPUT FILE %s",mtd_data->metaFilename);
        plumed_error(buf);
      }
    }
#else
    if(mtd_data->repl==-1){
      char buf[1024];
      sprintf(buf, "MISSING PLUMED INPUT FILE %s",mtd_data->metaFilename);
      plumed_error(buf);
    } else if(mtd_data->repl>-1) {
      strcpy(tmpmeta,mtd_data->metaFilename);
      tmpmeta[strlen(mtd_data->metaFilename) - 4] = '\0';
      sprintf(tmpmeta+strlen(tmpmeta),"%d",mtd_data->repl);
      sprintf(metafile, "%s.dat", tmpmeta);
      file = fopen(metafile, "r");
      if(file==NULL) {
        char buf[1024];
        sprintf(buf, "MISSING PLUMED INPUT FILE %s",mtd_data->metaFilename);
        plumed_error(buf);
      }
    }
#endif
  }

// CV counter initialization and calling routine to set the default values for many variables
  count = 0;
  iline = 0;
  read_defaults();

// Initialize everything for reconnaissance metadynamics
#ifdef RECONMETAD 
//  reconinpt.monitor=0; 
  reconinpt.nconst=0; reconinpt.nscales=0;
  reconinpt.stride=-10; reconinpt.height=-10.0; reconinpt.width=-10.0; 
  reconinpt.basTol=-10.0; reconinpt.diffconst=-10.0; reconinpt.iSize=-10.0;
  reconinpt.restart=0; reconinpt.fuzzy=0;
#endif 
// And initialize counter for bespoke collective coordinates
  colvar.nbespoke=0;  


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                               INPUT PARSER
//.............................................................................
// first word must be keyword (like COORD), then the parser
// seeks on the same line for additional input (like SIGMA).
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  fprintf(mtd_data->fplog, " \n");  
  fprintf(mtd_data->fplog, "::::::::::::::::: READING PLUMED INPUT :::::::::::::::::\n");

  plumed_read_input(&input,file,mtd_data->fplog);
// everything is now in the "input" structure, the actual file is not needed anymore
  fclose(file);

  for(iline=0;iline<input.nlines;iline++){

    nw   = input.nwords[iline];
    word = input.words[iline];

// empty line
    if(nw==0) continue;

// explicit comment, to be copied on the log file
    if(!strcmp(word[0],"NOTE") || !strcmp(word[0],"COMMENT")){
      fprintf(mtd_data->fplog, "\nCOMMENT: ");
      for(i=1;i<nw;i++)fprintf(mtd_data->fplog, "%s ", word[i]);
      fprintf(mtd_data->fplog,"\n");
// untested features
    } else if(!strcmp(word[0],"ENABLE_UNTESTED_FEATURES")){
      fprintf(mtd_data->fplog, "|-##################################################################\n");
      fprintf(mtd_data->fplog, "|- ENABLE_UNTESTED_FEATURES\n");
      fprintf(mtd_data->fplog, "|- THIS FLAG ENABLES FEATURES WHICH ARE NOT EXPLAINED IN THE MANUAL\n");
      fprintf(mtd_data->fplog, "|- AND COULD BE BUGGY\n");
      fprintf(mtd_data->fplog, "|- USE IT ONLY IF YOU ARE A PLUMED DEVELOPERS\n");
      fprintf(mtd_data->fplog, "|-##################################################################\n");
      logical.enable_untested_features=1;
    } else if(!strcmp(word[0],"NEW_COLVAR_FMT")){
      fprintf(mtd_data->fplog, "|- NEW COLVAR FMT ENABLED\n");
      mtd_data->newcolvarfmt=1;
    } else if(!strcmp(word[0],"OLD_COLVAR_FMT")){
      fprintf(mtd_data->fplog, "|- OLD COLVAR FMT ENABLED\n");
      mtd_data->newcolvarfmt=0;
// committors analysis
    } else if(!strcmp(word[0],"COMMITMENT")){
      logical.commit = 1;
      fprintf(mtd_data->fplog, "|-COMMITMENT ANALYSIS: YOU WILL ONLY MONITOR YOUR CVs MICRODYNAMICS\n");
      for(iw=1;iw<nw;iw++) if(!strcmp(word[iw],"NCV")){iw++; sscanf(word[iw],"%d",&commit.ncv);}
      for(iw=1;iw<nw;iw++){
       if(!strcmp(word[iw],"CV")){
          iw++; for(i=0;i<commit.ncv;i++) {sscanf(word[iw], "%d", &(commit.index[i])); commit.index[i] -= 1; iw++;};
          iw--;
       }     
      }
      for(i=0;i<commit.ncv;i++){
        iline++; // this mimicks a fgets
        sscanf(input.words[iline][0],"%lf",&uno);
        sscanf(input.words[iline][1],"%lf",&due);
        sscanf(input.words[iline][2],"%lf",&tre);
        sscanf(input.words[iline][3],"%lf",&quattro);
        ix=commit.index[i];
        commit.Amin[ix] = (real) uno;
        commit.Amax[ix] = (real) due;
        commit.Bmin[ix] = (real) tre;
        commit.Bmax[ix] = (real) quattro;
        fprintf(mtd_data->fplog, "|--CV %i: A min %f, max %f -- B min %f, max %f\n", ix+1, commit.Amin[ix], commit.Amax[ix], commit.Bmin[ix], commit.Bmax[ix]);
      }
      fprintf(mtd_data->fplog, "\n");
// Reconnaissance metadynamics
#ifdef RECONMETAD 
    } else if(!strcmp(word[0],"RECONNAISSANCE")){
      reconOn=1;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"RESTART")){ reconinpt.restart=1; logical.append=1; }
        else if(!strcmp(word[iw],"CV_LIST")){ iw++; reconinpt.nconst=plumed_get_group(word[iw],&reconinpt.cvlist,0,&input,mtd_data->fplog); }
//        else if(!strcmp(word[iw],"MONITOR")){ reconinpt.monitor=1; }
        else{ plumed_error("Unknown option for RECONNAISSANCE keyword"); }
      }
    } else if(!strcmp(word[0],"ONIONS")) {
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"HEIGHT")){
           iw++; sscanf(word[iw],"%lf",&uno); reconinpt.height = (real) uno*mtd_data->eunit;
        } else if(!strcmp(word[iw],"W_STRIDE")){
          iw++; sscanf(word[iw],"%i",&reconinpt.stride);
        } else if(!strcmp(word[iw],"WIDTH")){
          iw++; sscanf(word[iw],"%lf",&uno); reconinpt.width = (real) uno;
        } else {
          plumed_error("Unknown option for ONIONS keyword");
        }
      }
    } else if(!strcmp(word[0],"BASINS")) {
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"BASIN_TOL")){
           iw++; sscanf(word[iw],"%lf",&uno); reconinpt.basTol = (real) uno;
        } else if(!strcmp(word[iw],"EXPAND_PARAM")){
          iw++; sscanf(word[iw],"%lf",&uno); reconinpt.diffconst = (real) uno;
        } else if(!strcmp(word[iw],"INITIAL_SIZE")){
          iw++; sscanf(word[iw],"%lf",&uno); reconinpt.iSize = (real) uno;
        } else {
          plumed_error("Unknown option for BASINS keyword");
        }
      }
    } else if(!strcmp(word[0],"CLUSTER")) {
      reconinpt.runFreq[reconinpt.nscales]=-10; reconinpt.storeFreq[reconinpt.nscales]=-10;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"RUN_FREQ")){
           iw++; sscanf(word[iw],"%i",&(reconinpt.runFreq[reconinpt.nscales]) );
        } else if(!strcmp(word[iw],"STORE_FREQ")){
          iw++; sscanf(word[iw],"%i",&(reconinpt.storeFreq[reconinpt.nscales]) );
        } else if(!word[iw],"FUZZY"){
          reconinpt.fuzzy=1;
        } else {
          plumed_error("Unknown option for CLUSTER keyword");
        }
      }
      reconinpt.nscales++;
      //fprintf(mtd_data->fplog,"NUMBER OF CLUSTER SCALES %i \n",reconinpt.nscales);
#else
    } else if(!strcmp(word[0],"RECONNAISSANCE")){
      plumed_error("To run reconMetaD you must patch with lapack and lstdc++ libraries.  Revert and repatch");
    } else if(!strcmp(word[0],"ONIONS")) {
      plumed_error("To run reconMetaD you must patch with lapack and lstdc++ libraries.  Revert and repatch");
    } else if(!strcmp(word[0],"BASINS")) {
      plumed_error("To run reconMetaD you must patch with lapack and lstdc++ libraries.  Revert and repatch");
    } else if(!strcmp(word[0],"CLUSTER")) {
      plumed_error("To run reconMetaD you must patch with lapack and lstdc++ libraries.  Revert and repatch");
#endif
// metadynamics
    } else if(!strcmp(word[0],"HILLS")) {
      logical.do_hills = 1;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"HEIGHT")){
          iw++; sscanf(word[iw],"%lf",&uno); hills.wwr = (real) uno*mtd_data->eunit;
        } else if(!strcmp(word[iw],"W_STRIDE")){
          iw++; sscanf(word[iw],"%i",&hills.nt_hills);
        } else if(!strcmp(word[iw],"R_STRIDE")){
          iw++; sscanf(word[iw],"%i",&hills.nr_hills);
        } else if(!strcmp(word[iw],"RESTART")){
          logical.restart_hills = 1;
        } else if(!strcmp(word[iw],"RATE")){
          iw++; sscanf(word[iw],"%lf",&uno); hills.rate = (real) uno;
        } else if(!strcmp(word[iw],"MAX_HEIGHT")){
          iw++; sscanf(word[iw],"%lf",&uno); hills.max_height = (real) uno;
        } else if(!strcmp(word[iw],"MAX_STRIDE")){
          iw++; sscanf(word[iw],"%i",&hills.max_stride);
        } else {
          plumed_error("Unknown option for HILLS keyword");
        }
      };
      if(hills.rate==0.0) hills.rate = hills.wwr/hills.nt_hills/mtd_data->dt;
      if(hills.wwr==0.0) hills.wwr=hills.rate*hills.nt_hills*mtd_data->dt;
      fprintf(mtd_data->fplog,"|-HILLS:\n");
      if(hills.max_height==0.0){
        fprintf(mtd_data->fplog,"|--HEIGHT %f  WRITING STRIDE %i DEPOSITION RATE %f \n",
                hills.wwr/mtd_data->eunit, hills.nt_hills, hills.rate/mtd_data->eunit);
      } else {
        fprintf(mtd_data->fplog,"|--DEPOSITION RATE %f \n",hills.rate/mtd_data->eunit);
        fprintf(mtd_data->fplog,"|--MAXIMUM STRIDE BETWEEN HILLS %i \n",hills.max_stride);
        fprintf(mtd_data->fplog,"|--MAXIMUM HEIGHT               %f \n",hills.max_height);
      }
      if(logical.restart_hills) fprintf(mtd_data->fplog,"|-RESTARTING METADYNAMICS!\n");
      if(hills.nr_hills!=1) fprintf(mtd_data->fplog,"|--READING STRIDE %i\n", hills.nr_hills);
      fprintf(mtd_data->fplog, "\n");
// multiple walkers (shared file)
    } else if(!strcmp(word[0],"MULTIPLE_WALKERS")) {
      logical.do_walkers = 1; 
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"HILLS_DIR")){
          iw++; sscanf(word[iw], "%s", hills.dir);
        } else if(!strcmp(word[iw],"NWALKERS")){
          iw++; sscanf(word[iw], "%i", &hills.nwalkers); 
        } else if(!strcmp(word[iw],"ID")){
          iw++; sscanf(word[iw], "%i", &hills.idwalker);
        } else if(!strcmp(word[iw],"R_STRIDE")){
          iw++; sscanf(word[iw], "%i", &hills.nr_hills); 
        } else {
          plumed_error("Unknown option for MULTIPLE_WALKERS keyword");
        }
      };
      fprintf(mtd_data->fplog,"|-MULTIPLE WALKERS: NWALKERS %i ID %i\n", hills.nwalkers,hills.idwalker);
      if(hills.nr_hills!=1) fprintf(mtd_data->fplog,"|--READING STRIDE %i\n", hills.nr_hills);
      fprintf(mtd_data->fplog,"|--DIRECTORY FOR HILLS I/O %s\n", hills.dir);
    } else if(!strcmp(word[0],"PRINT")) {
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"W_STRIDE")){
          iw++; sscanf(word[iw],"%i", &colvar.nt_print);
        } else if(!strcmp(word[iw],"T_OFFSET")){
          iw++; sscanf(word[iw],"%lf", &uno); mtd_data->time_offset= uno;  
        } else if(!strcmp(word[iw],"APPEND")){
          logical.append = 1;
        } else {
          plumed_error("Unknown option for PRINT keyword");
        }
      }
      fprintf(mtd_data->fplog,"|-PRINTING ON COLVAR FILE EVERY %i STEPS\n",colvar.nt_print);
      fprintf(mtd_data->fplog,"|-INITIAL TIME OFFSET IS %f TIME UNITS\n",mtd_data->time_offset);
      if(logical.append) fprintf(mtd_data->fplog,"|-APPENDING TO OLD COLVAR\n");
      logical.print = 1;
    } else if(!strcmp(word[0],"HILLS_LABEL")) {
      if( nw==1 ){
        plumed_error("Missing label after HILLS_LABEL"); 
      }else{
        sscanf(word[1],"%s",colvar.hills_label);
      }
      fprintf(mtd_data->fplog,"|-USING HILLS_LABEL %s IN HEADERS OF COLVAR AND HILLS FILES\n",colvar.hills_label);
    } else if(!strcmp(word[0],"DUMP_ATOMS")) {
#if ! defined (PLUMED_GROMACS)
        plumed_error("DUMP_ATOMS NOT YET IMPLEMENTED IN THIS CODE");
#endif
       int list_found;
       list_found=0;
       fprintf(mtd_data->fplog,"|- DUMPING ATOMS\n");
       for(iw=1;iw<nw;iw++){
         if(!strcmp(word[iw],"LIST")){
           list_found=1;
           iw++; mtd_data->dump_atoms+=plumed_get_group(word[iw],&mtd_data->dump_list,mtd_data->dump_atoms,&input,mtd_data->fplog);
         } else if(!strcmp(word[iw],"STRIDE")){
           iw++; sscanf(word[iw],"%i", &mtd_data->dump_stride);
         } else {
           plumed_error("Unknown option for DUMP keyword");
         }
       }
       if(!list_found)plumed_error("NEEDED LIST KEYWORD FOR DUMP_ATOMS\n");
       fprintf(mtd_data->fplog,"|- SET MEMBERS: ");
       for(i=0;i<mtd_data->dump_atoms;i++){
         fprintf(mtd_data->fplog," %d ",mtd_data->dump_list[i]+1);if((i+1)%20==0)fprintf(mtd_data->fplog,"\n               ");
       }fprintf(mtd_data->fplog,"\n\n");
     } else if(!strcmp(word[0],"ALIGN_ATOMS")) {
       int list_found;
       list_found=0;
       fprintf(mtd_data->fplog,"|- ALIGNING ATOMS\n");
       for(iw=1;iw<nw;iw++){
         if(!strcmp(word[iw],"LIST")){
           list_found=1;
           iw++; colvar.align_atoms+=plumed_get_group(word[iw],&colvar.align_list,colvar.align_atoms,&input,mtd_data->fplog);
         } else {
           plumed_error("Unknown option for ALIGN keyword");
         }
       }
       if(!list_found)plumed_error("NEEDED LIST KEYWORD FOR ALIGN_ATOMS\n");
       fprintf(mtd_data->fplog,"|- SET MEMBERS: ");
       for(i=0;i<colvar.align_atoms;i++){
         fprintf(mtd_data->fplog," %d ",colvar.align_list[i]+1);if((i+1)%20==0)fprintf(mtd_data->fplog,"\n               ");
       }fprintf(mtd_data->fplog,"\n\n");
    } else if(!strcmp(word[0],"DEBUG_DERIVATIVES")){
      logical.debug_derivatives=1;
    } else if(!strcmp(word[0],"PARALLEL_HILLS")){
#if ! defined (PLUMED_GROMACS) && ! defined (DL_POLY) && ! defined (AMBER)
        plumed_error("PARALLEL_HILLS NOT YET IMPLEMENTED IN THIS CODE");
#endif
        for(iw=1;iw<nw;iw++){
          if(!strcmp(word[iw],"ON")){
            logical.parallel_hills=1;
          } else if(!strcmp(word[iw],"OFF")){
            logical.parallel_hills=0;
          } else {
            plumed_error("Unknown flag for keyword PARALLEL_HILLS");
          }
        }
        fprintf(mtd_data->fplog, "|-PARALLEL HILLS ");
        if(logical.parallel_hills) fprintf(mtd_data->fplog, "ON");
        else                       fprintf(mtd_data->fplog, "OFF");
        fprintf(mtd_data->fplog, "\n\n");
    } else if(!strcmp(word[0],"PTMETAD")){
      #if ! defined (PLUMED_GROMACS) && ! defined (OPEP)
          plumed_error("PTMETAD: NOT YET IMPLEMENTED IN THIS CODE");
      #else
      fprintf(mtd_data->fplog, "|-PARALLEL TEMPERING METADYNAMICS\n");
      logical.remd = 1;
      fprintf(mtd_data->fplog, "|--REPLICA 0 TEMPERATURE = %f\n", mtd_data->rte0);
      fprintf(mtd_data->fplog, "|--REPLICA %i TEMPERATURE = %f\n", mtd_data->repl, mtd_data->rteio);
      if(mtd_data->repl==-1) {
        fprintf(mtd_data->fplog, "\n!!!! mdrun not in replica exchange mode, keyword PTMETAD will not be considered !!!!\n");
        logical.remd = 0;
      }
      fprintf(mtd_data->fplog, "\n");
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"NEIGHBOUR")){
          iw++; sscanf(word[iw],"%i",&colvar.ptmetad_neighbours);
        } else if(!strcmp(word[iw],"SIGMA")){
          iw++; sscanf(word[iw],"%lf",&uno); colvar.ptmetad_sigma = (real) uno;
        } else if(!strcmp(word[iw],"HREX")){
#if ! defined (PLUMED_GROMACS)
          plumed_error("HAMILTONIAN REPLICA-EXCHANGE: NOT YET IMPLEMENTED IN THIS CODE");
#endif
          fprintf(mtd_data->fplog, "|-HAMILTONIAN REPLICA-EXCHANGE\n");
          logical.hrex=1;
        } else if(!strcmp(word[iw],"NORESCALE")){
          logical.norescale=1;
          fprintf(mtd_data->fplog, "|--DO NOT RESCALE GAUSSIANS HEIGHT WITH TEMPERATURE\n"); 
        } else {
          plumed_error("Unknown flag for keyword PTMETAD");
        }
      }
      if(colvar.ptmetad_neighbours){
        fprintf(mtd_data->fplog, "|--SIGMA = %lf\n",(double) colvar.ptmetad_sigma);
        fprintf(mtd_data->fplog, "|--NEIGHBOUR = %i\n",colvar.ptmetad_neighbours);
      }
      fprintf(mtd_data->fplog, "\n");
      #endif
    } else if(!strcmp(word[0],"BIASXMD")){
      #if ! defined (PLUMED_GROMACS) && ! defined (OPEP)
          plumed_error("|-BIASXMD: NOT YET IMPLEMENTED IN THIS CODE");
      #else
      fprintf(mtd_data->fplog, "|-BIAS EXCHANGE METADYNAMICS\n");
      logical.rpxm = 1;
      logical.remd = 1;
      if(mtd_data->repl==-1&&!(mtd_data->mcr->ms)) {
        plumed_error("\n!!!! mdrun not in replica exchange mode, keyword BIASXMD cannot be used !!!!\n");
        logical.rpxm = 0;
        logical.remd = 0;
      }
      fprintf(mtd_data->fplog, "\n");
      #endif
    } else if(!strcmp(word[0],"TAMD") ){
      // ### Removed because keyword DAFED is reserved.
      //  || !strcmp(word[0],"DAFED")
      logical.tamd = 1;
      int read_biasfactor = 0;
      int read_cvtemp = 0;
      int read_simtemp = 0;
      int read_starttemp = 0;
      tamd.seed=1234;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CVTEMP")){
          iw++; sscanf(word[iw], "%lf", &uno); tamd.wtemp   = (real) uno; read_cvtemp = 1;
        } else if(!strcmp(word[iw],"TFACTOR")){
          iw++; sscanf(word[iw], "%lf", &uno); tamd.wfactor = (real) uno; read_biasfactor = 1;
        } else if(!strcmp(word[iw],"SIMTEMP")){
          iw++; sscanf(word[iw], "%lf", &uno); tamd.simtemp = (real) uno; read_simtemp = 1;
        } else if(!strcmp(word[iw],"STARTTEMP")){
          iw++; sscanf(word[iw], "%lf", &uno); tamd.starttemp = (real) uno; read_starttemp = 1;
        } else if(!strcmp(word[iw],"TAU")){
          iw++; sscanf(word[iw], "%lf", &uno); tamd.tau= (real) uno;
        } else if(!strcmp(word[iw],"SEED")){
          iw++; sscanf(word[iw], "%i", &tamd.seed); 
        } else {
          plumed_error("Unknown flag for keyword TAMD/DAFED");
        }
      }
      tamd.drift=0.0;
      if(!read_simtemp)
        plumed_error("WITH TAMD/DAFED YOU ALWAYS HAVE TO SPECIFY THE \"SIMTEMP \" KEYWORD");
      if(read_biasfactor==read_cvtemp)
        plumed_error("WITH TAMD/DAFED YOU HAVE TO SPECIFY EITHER \"CVTEMP \" OR \"TFACTOR \" KEYWORD");
      if(read_cvtemp)     tamd.wfactor = tamd.wtemp / tamd.simtemp;  
      if(read_biasfactor) tamd.wtemp   = tamd.wfactor * tamd.simtemp;
      if(!read_starttemp) tamd.starttemp = tamd.simtemp;
      fprintf(mtd_data->fplog, "|-TAMD/DAFED WITH TFACTOR %f (CVTEMP = %f) \n", tamd.wfactor,tamd.wtemp);
      fprintf(mtd_data->fplog, "|-SPRING CONSTANTS CALCULATED FROM SIGMA AND WTEMP OF EACH VARIABLE\n");
      fprintf(mtd_data->fplog, "|-INITIAL RANDOMIZATION OF RESTRAINTS AT T = %f\n",tamd.simtemp);
      fprintf(mtd_data->fplog, "\n");
    } else if(!strcmp(word[0],"WELLTEMPERED")){
      logical.welltemp = 1;
      int read_biasfactor = 0;
      int read_cvtemp = 0;
      int read_simtemp = 0;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CVTEMP")){
          iw++; sscanf(word[iw], "%lf", &uno); colvar.wtemp   = (real) uno; read_cvtemp = 1;
        } else if(!strcmp(word[iw],"BIASFACTOR")){
          iw++; sscanf(word[iw], "%lf", &uno); colvar.wfactor = (real) uno; read_biasfactor = 1;
        } else if(!strcmp(word[iw],"SIMTEMP")){
          iw++; sscanf(word[iw], "%lf", &uno); colvar.simtemp = (real) uno; read_simtemp = 1;
        } else if(!strcmp(word[iw],"READ_OLD_BF")){
          logical.read_old_bf = 1;
        // JFD>
        } else if(!strcmp(word[iw],"DICKSONIAN_TEMPERING")){
          logical.dicksonian_tempering = 1;
        // <JFD
        } else {
          plumed_error("Unknown flag for keyword WELLTEMPERED");
        }
      }
      if(!read_simtemp)
        plumed_error("WITH WELLTEMPERED YOU ALWAYS HAVE TO SPECIFY THE \"SIMTEMP \" KEYWORD");
      if(read_biasfactor==read_cvtemp)
        plumed_error("WITH WELLTEMPERED YOU HAVE TO SPECIFY EITHER \"CVTEMP \" OR \"BIASFACTOR \" KEYWORD");
      if(read_cvtemp)     colvar.wfactor = colvar.wtemp / colvar.simtemp;  
      if(read_biasfactor) colvar.wtemp   = colvar.wfactor * colvar.simtemp;
      if (colvar.wfactor<=1.0) {  
        char buf[1024];
        sprintf(buf,"WELLTEMPERED: BIASFACTOR less than or equal to 1.0 ( %f ) \n",colvar.wfactor);
        plumed_error(buf);
      } 
      fprintf(mtd_data->fplog, "|-WELL-TEMPERED METADYNAMICS WITH BIASFACTOR %f (CVTEMP = %f) \n", colvar.wfactor,colvar.wtemp);
      if(logical.read_old_bf) fprintf(mtd_data->fplog, "|--READING OLD BIASFACTOR WHEN RESTARTING \n");
      if(logical.dicksonian_tempering) fprintf(mtd_data->fplog, "|--USING DICKSONIAN TEMPERING\n");
      fprintf(mtd_data->fplog, "\n"); 
    // JFD>
    } else if(!strcmp(word[0],"TRANSITIONTEMPERED")){
      logical.transition_tempering = 1;
      int read_biasfactor = 0;
      int read_cvtemp = 0;
      int read_simtemp = 0;
      int ncv;
      int read_ncv = 0;
      int read_wells = 0;
      colvar.ttthreshold = 0.0;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CVTEMP")){
          iw++; sscanf(word[iw], "%lf", &uno); colvar.tttemp   = (real) uno; read_cvtemp = 1;
        } else if(!strcmp(word[iw],"BIASFACTOR")){
          iw++; sscanf(word[iw], "%lf", &uno); colvar.ttfactor = (real) uno; read_biasfactor = 1;
        } else if(!strcmp(word[iw],"SIMTEMP")){
          iw++; sscanf(word[iw], "%lf", &uno); colvar.simtemp = (real) uno; read_simtemp = 1;
        } else if(!strcmp(word[iw],"NCV")){
          iw++; sscanf(word[iw], "%d", &ncv); read_ncv = 1;
          if (ncv > nconst_max) plumed_error("Invalid choice for NCV for keyword TRANSITIONTEMPERED");
        } else if(!strcmp(word[iw],"TRANSITIONWELL")){
          if (!read_ncv) plumed_error("Must specify NCV before TRANSITIONWELL for keyword TRANSITIONTEMPERED");
          for (i=0;i<ncv;i++) {
            iw++; sscanf(word[iw], "%lf", &uno);
            colvar.transition_wells[read_wells][i] = (real) uno;
          }
          read_wells++;
        } else if(!strcmp(word[iw],"THRESHOLD")){
          iw++; sscanf(word[iw], "%lf", &uno);
          colvar.ttthreshold = mtd_data->boltz * (real) uno;
        } else if(!strcmp(word[iw],"READ_OLD_BF")){
          logical.read_old_bf = 1;
        } else {
          plumed_error("Unknown flag for keyword TRANSITIONTEMPERED");
        }
      }
      if(!(read_wells >= 2))
        plumed_error("WITH TRANSITIONTEMPERED YOU ALWAYS HAVE TO SPECIFY AT LEAST TWO WELLS WITH THE \"TRANSITIONWELL \" KEYWORD");
      colvar.n_transition_wells = read_wells;
      if(!read_simtemp)
        plumed_error("WITH TRANSITIONTEMPERED YOU ALWAYS HAVE TO SPECIFY THE \"SIMTEMP \" KEYWORD");
      if(read_biasfactor==read_cvtemp)
        plumed_error("WITH TRANSITIONTEMPERED YOU HAVE TO SPECIFY EITHER \"CVTEMP \" OR \"BIASFACTOR \" KEYWORD");
      if(read_cvtemp)     colvar.ttfactor = colvar.tttemp / colvar.simtemp;  
      if(read_biasfactor) colvar.tttemp   = colvar.ttfactor * colvar.simtemp;
      if (colvar.ttfactor<=0.0) {  
        char buf[1024];
        sprintf(buf,"TRANSITIONTEMPERED: BIASFACTOR less than or equal to 0.0 ( %f ) \n",colvar.wfactor);
        plumed_error(buf);
      } 
      fprintf(mtd_data->fplog, "|-TRANSITION-TEMPERED METADYNAMICS WITH BIASFACTOR %f (CVTEMP = %f) \n", colvar.wfactor,colvar.wtemp);
      if(logical.read_old_bf) fprintf(mtd_data->fplog, "|--READING OLD BIASFACTOR WHEN RESTARTING \n");
      fprintf(mtd_data->fplog, "\n");
  // <JFD
  // ADW>
    }else if(!strcmp(word[0],"TARGET_DISTRIBUTION")) {
      logical.target_distribution = 1;
      fprintf(mtd_data->fplog,"Enabling target distribution metadynamics\n");
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"FILENAME")){
          iw++; 
	  sscanf(word[iw], "%s", target_grid.r_file);
        }
        if(!strcmp(word[iw],"SIMTEMP")){
          iw++; 
	  sscanf(word[iw], "%lf", &uno);
	  colvar.simtemp = uno;
        }
      }
    } else if(!strcmp(word[0], "EDS")) {
      //EXAMPLE for adaptive
      //EDS STRIDE 500 SIMTEMP 300 SEED 4143 FILENAME FOO CV LIST 1 2 4
      //EDS CV CENTERS 0.5 2.5 2.3
      //EDS CV RANGES 5 5 5
      //EXAMPLE for fixed
      //EDS SIMTEMP 300 CV LIST 1 2 4
      //EDS CV CENTERS 0.5 2.5 2.3
      //EDS CV RANGES 34 23 4      
      if(!logical.eds) 
	fprintf(mtd_data->fplog, "Enabling experiment directed simulation\n");
      logical.eds = 1;
      iw = 1;
      
      if(!strcmp(word[iw], "CV")) {
	iw++;
	if(!strcmp(word[iw++], "CENTERS")) {
	  if(eds.cv_number == 0) {
	    plumed_error("Must define CVs first for EDS with [CV LIST 1 2 3]");
	  }
	  for(icv = 0;iw < nw; iw++) {
	    sscanf(word[iw], "%lf", &uno);
	    eds.centers[icv] = uno;
	    fprintf(mtd_data->fplog, "EDS: Will center CV %d at %lf\n", icv+1,
		    eds.centers[icv]);
	    icv++;
	  }
	}
	else if(!strcmp(word[iw - 1], "RANGES")) {
	  if(eds.cv_number == 0) {
	    plumed_error("Must define CVs first for EDS with [CV LIST 1 2 3]");
	  }
	  for(icv = 0;iw < nw; iw++) {
	    sscanf(word[iw], "%lf", &uno);
	    eds.max_coupling_range[icv] = uno;
	    eds.max_coupling_rate[icv] = eds.max_coupling_range[icv] / (10 * eds.update_period);
	    fprintf(mtd_data->fplog, 
		    "EDS: Will cap range of CV %d at %lf\n", 
		    icv+1,
		    eds.max_coupling_range[icv]);
	    icv++;
	  }
	} else {
	  plumed_error("Syntax is EDS CV RANGES.... or EDS CV CENTERS....\n");
	}
      } else {

	if(eds.cv_number != 0) {
	  plumed_error("Syntax is EDS CV RANGES.... or EDS CV CENTERS....\n");
	}

	if(iw >= nw - 2 || strcmp(word[iw++], "STRIDE"))
	  plumed_error("Must specify STRIDE in EDS [EDS 500 CV LIST 1 3]\n");
	int update_period;
	if(!sscanf(word[iw++], "%d", &update_period)){
	  plumed_error("Must specify STRIDE in EDS [EDS STRIDE 500 SIMTEMP 300 SEED 431 CV LIST 1 3]\n");

	}

	if(iw >= nw - 2 || strcmp(word[iw++], "SIMTEMP"))
	  plumed_error("Must specify SIMTEMP in EDS [EDS STRIDE 500 SIMTEMP 300 SEED 431 CV LIST 1 3]\n");
	if(!sscanf(word[iw++], "%lf", &uno)){
	  plumed_error("Must specify SIMTEMP in EDS [EDS STRIDE 500 SIMTEMP 300 SEED 4313 CV LIST 1 3]\n");
	}

	int eds_seed = 0;
	if(!strcmp(word[iw], "SEED")) {
	  if(!sscanf(word[++iw], "%d", &eds_seed)){
	    plumed_error("Must use integer SEED\n");
	  }
	  iw++;
	}

	char* filename = "EDS_OUT";
	if(!strcmp(word[iw], "FILENAME")) {
	  if(!sscanf(word[iw++], "%s", filename)){
	    plumed_error("Filename invalid\n");
	  }
	  iw++;
	}
		
	if(iw < nw - 3 || !(strcmp(word[iw++], "CV") & strcmp(word[iw++], "LIST"))) {
	  int* cv_map = (int*) malloc(sizeof(int) * colvar.nconst);
	  for(i = 0;iw < nw; iw++) {
	    sscanf(word[iw], "%d", &icv);
	    cv_map[i] = icv - 1;
	    fprintf(mtd_data->fplog, 
		    "EDS: Will use CV %d \n", 
		    icv);
	    i++;
	  }
	  cv_map = (int *) realloc(cv_map, sizeof(int) * i);
	  eds_init(i, update_period, uno, eds_seed, 0, cv_map, (const char*) filename, &eds);
	} else {
	  plumed_error("Must specify CV List in EDS [EDS 500 CV LIST 1 3]\n");
	}
      }
      
    }else if(!strcmp(word[0],"TREAT_INDEPENDENT")) {
      //Easy to test here
      if(!logical.enable_untested_features)
	plumed_error("TREAT_INDEPENDENT NOT ENABLED!\n");

      //parse list of CVs
      fprintf(mtd_data->fplog,"Will sample CVs ");
      for(iw = 1; iw < nw; iw++) {
	//we assume it's a stochastic sample setting first, then
	//check if it's actually an integer
	if(sscanf(word[iw], "%lf", &uno)) {
	  if(ceil(uno) == uno) {
	    icv = (int) uno;
	    colvar.b_treat_independent[icv-1] = 1;
	    fprintf(mtd_data->fplog,"%d ", icv);
	  } else {
	    break;
	  }
	} else {
	  plumed_error("INCORRECT FORMAT FOR INDEPENDENT. SHOULD BE: \n\
                        TREAT_INDEPENDENT (CV IDS) [STOCHASTIC SAMPLING FRACTION] [SEED FOR SAMPLING]\
                        \nEXAMPLE:\n TREAT_INDEPENDENT 1 2 3 0.25 43323");
	}
      }
      fprintf(mtd_data->fplog," independently\n");
      
      //Is the sampling going to be stochastic?
      if(ceil(uno) != uno) {
	if(uno >= 1 || uno < 0) 
	  plumed_error("SAMPLING FRACTION MUST BE BETWEEN 0 AND 1\n");
	colvar.stoch_sample = uno;	     
 	fprintf(mtd_data->fplog,"Will sample independent CVs with probability %lf \n", colvar.stoch_sample);
      }
      
      //Is the seed specified?
      if(iw >= nw || !sscanf(word[iw], "%d", &colvar.stoch_sample_seed)) {
	colvar.stoch_sample_seed = 0;
	fprintf(mtd_data->fplog, "DEFAULT STOCHASTIC SAMPLING SEED IS 0\n");
      }	else
	fprintf(mtd_data->fplog,"Will use seed %d for sampling\n", colvar.stoch_sample_seed);
      
      // <ADW
    } else if(!strcmp(word[0],"DEBUG_TRANSITIONTEMPERED")){
      logical.ttdebug=1;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"FILENAME")){
          iw++; sscanf(word[iw], "%s", colvar.ttdebug_file);
        }
      };
      fprintf(mtd_data->fplog, "|-WRITING DEBUG_TRANSITIONTEMPERED ON FILE \n");
      fprintf(mtd_data->fplog, "|--FILENAME %s\n\n", colvar.ttdebug_file);
    } else if(!strcmp(word[0],"DEBUG")){
      fprintf(mtd_data->fplog,"|- CV DERIVATIVES DEBUGGING \n"); 
      logical.debug = 1;
    } else if(!strcmp(word[0],"DISTANCE")){
      read_dist(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"MINDIST")){
      read_mindist(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"COORD")){
      read_coord(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"ANGLE")){
      read_angle(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"HBONDS")){
      read_hbonds(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"TORSION")){
      read_torsion(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"RGYR")){
      read_rgyr(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"DIPOLE")){
      read_dipole(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"DIHCOR")) {
      read_dihcor(word, count, &input,&iline,mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"WATERBRIDGE")) {
      read_waterbridge(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"DENSITY")) {
      read_density(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"DENSITYSWITCH")) {
      read_densityswitch(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"ALPHABETA")) {
      read_alfabeta(word, count,&input,&iline,mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"S_PATH")) {
      logical.path = 1;
      colvar.type_s[count]   = 30;
      read_path(word, count, &input,mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"Z_PATH")) {
      logical.path = 1;
      colvar.type_s[count]   = 31;
      read_path(word, count, &input,mtd_data->fplog);
      count++;   
    } else if(!strcmp(word[0],"TARGETED")) {
      colvar.type_s[count]   = 31;
      read_path(word, count, &input,mtd_data->fplog);
      count++;   
/* Atom position */      
    } else if(!strcmp(word[0],"POSITION")) {
      colvar.type_s[count]   = 32;
      read_position(word, count, &input, mtd_data->fplog);
      count++;   
    } else if(!strcmp(word[0],"ELSTPOT")) {
      colvar.type_s[count]   = 33;
      read_elstpot(word, count, &input, mtd_data->fplog);
      count++;   
    } else if(!strcmp(word[0],"PUCKERING")) {
      logical.puckering      =  1;
      colvar.type_s[count]   = 34;
      read_puckering(word, count, &input, mtd_data->fplog);
      count++;   
    } else if(!strcmp(word[0],"ENERGY")) {
#if ! defined (PLUMED_GROMACS4) && ! defined (DL_POLY) && ! defined (AMBER) && ! defined (PLUMED_GROMACS45)
          plumed_error("ENERGY CV: NOT YET IMPLEMENTED IN THIS CODE");
      #endif
      logical.energy         =  1;
      colvar.type_s[count]   = 35;
      read_energy(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"HELIX")) {
      colvar.type_s[count]   = 36;
      read_helix(word, count,&input,&iline,mtd_data->fplog);
      count++; 
    } else if(!strcmp(word[0],"ALPHARMSD")) {
      colvar.type_s[count]   = 37;
      read_alpharmsd(word, count,&input, mtd_data->fplog);
      logical.do_alphabetarmsd=1;
      count++;
    } else if(!strcmp(word[0],"ANTIBETARMSD")) {
      colvar.type_s[count]   = 38;
      read_antibetarmsd(word, count,&input, mtd_data->fplog);
      logical.do_alphabetarmsd=1;
      count++;
    } else if(!strcmp(word[0],"PARABETARMSD")) {
      colvar.type_s[count]   = 39;
      read_parabetarmsd(word, count,&input, mtd_data->fplog);
      logical.do_alphabetarmsd=1;
      count++;
    /*} else if(!strcmp(word[0],"CAMSHIFT")) {
      colvar.type_s[count]   = 40;
      read_camshift(word, count,&input, mtd_data->fplog);
      count++;*/
    } else if (!strcmp(word[0],"PCA")) {
      logical.do_pca = 1;
      colvar.type_s[count]   = 42;
      read_pca(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"CMAP")) {
      colvar.type_s[count]   = 45;
      read_cmap(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"POLY")) {
      colvar.type_s[count]   = 50;
      read_poly(word, count, &input, &iline, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"FUNCTION")) {
      colvar.type_s[count]   = 51;
      read_func(word, count, &input, &iline, &nw,  mtd_data->fplog);
      count++;
#ifdef CVS
    } else if(!strcmp(word[0],"BESPOKE")) {
      colvar.type_s[count] = 46;
      read_bespoke(word, count, &input, mtd_data->fplog);
      count++; colvar.nbespoke++;
#else   
    } else if(!strcmp(word[0],"BESPOKE")) {
      plumed_error("To run bespoke collective coordinates you must patch with lapack and lstc++ libraries.  Revert and repatch");
#endif
    } else if(!strcmp(word[0],"RDF")) {
      colvar.type_s[count] = 47;
      read_rdf(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"HISTOGRAM")) {
      colvar.type_s[count] = 49;
      read_histogram(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"ADF")) {
      colvar.type_s[count] = 52;
      read_adf(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"SPRINT")) {
      colvar.type_s[count]   = 55;
      read_sprint(word, count,&input, mtd_data->fplog);
      logical.do_sprint=1;
      count++;
    } else if(!strcmp(word[0],"MSD")) {
		// this is a simple dummy MSD variable whose intent is to be used into 
		// hybrid path and hypothetical should contain a more advanced path structure
		colvar.type_s[count] = 53;
		read_msd(word, count, &input, mtd_data->fplog);
		count++;
    } 
	  else if(!strcmp(word[0],"UWALL")){
      int read_limit=0;
      int read_kappa=0;
// first we select the proper CV
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv);}
      else{plumed_error("WITH UWALL YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");}
// then we parse the line
      logical.upper[icv-1]=1;
      logical.do_walls = 1;   // ### For modified output format
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")) {
          iw++;  // already read
        } else if(!strcmp(word[iw],"LIMIT")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.upper[icv-1]=(real)uno; read_limit=1;
        } else if(!strcmp(word[iw],"KAPPA")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.sigma[icv-1]=(real)uno*mtd_data->eunit; read_kappa=1;
        } else if(!strcmp(word[iw],"EXP")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.uexp[icv-1]=(real)uno;
        } else if(!strcmp(word[iw],"EPS")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.ueps[icv-1]=(real)uno;
        } else if(!strcmp(word[iw],"OFF")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.uoff[icv-1]=(real)uno;
        } else {
          plumed_error("Unknown flag for keyword UWALL");
        };
      }
      if(!read_limit)
        plumed_error("WITH UWALL YOU ALWAYS HAVE TO SPECIFY THE \"LIMIT\" KEYWORD\n");
      if(!read_kappa)
        plumed_error("WITH UWALL YOU ALWAYS HAVE TO SPECIFY THE \"KAPPA\" KEYWORD\n");
      fprintf(mtd_data->fplog, "|-WALL ON COLVAR %i: UPPER LIMIT = %f, KAPPA = %f, EXPONENT = %i, REDUX = %f, OFFSET = %f \n\n",
             icv, cvw.upper[icv-1], cvw.sigma[icv-1]/mtd_data->eunit, cvw.uexp[icv-1], cvw.ueps[icv-1], cvw.uoff[icv-1]);
     } else if(!strcmp(word[0],"LWALL")){
      int read_limit=0;
      int read_kappa=0;
// first we select the proper CV
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv);}
      else{plumed_error("WITH UWALL YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");}
// then we parse the line
      logical.lower[icv-1]=1;
      logical.do_walls = 1;   // ### For modified output format
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")) {
          iw++;  // already read
        } else if(!strcmp(word[iw],"LIMIT")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.lower[icv-1]=(real)uno; read_limit=1;
        } else if(!strcmp(word[iw],"KAPPA")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.lsigma[icv-1]=(real)uno*mtd_data->eunit; read_kappa=1;
        } else if(!strcmp(word[iw],"EXP")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.lexp[icv-1]=(real)uno;
        } else if(!strcmp(word[iw],"EPS")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.leps[icv-1]=(real)uno;
        } else if(!strcmp(word[iw],"OFF")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.loff[icv-1]=(real)uno;
        } else {
          plumed_error("Unknown flag for keyword LWALL");
        };
      }
      if(!read_limit)
        plumed_error("WITH LWALL YOU ALWAYS HAVE TO SPECIFY THE \"LIMIT\" KEYWORD\n");
      if(!read_kappa)
        plumed_error("WITH LWALL YOU ALWAYS HAVE TO SPECIFY THE \"KAPPA\" KEYWORD\n");
      fprintf(mtd_data->fplog, "|-WALL ON COLVAR %i: LOWER LIMIT = %f, KAPPA = %f, EXPONENT = %i, REDUX = %f, OFFSET = %f \n\n",
             icv, cvw.lower[icv-1], cvw.lsigma[icv-1]/mtd_data->eunit, cvw.lexp[icv-1], cvw.leps[icv-1], cvw.loff[icv-1]);
    // JFD> 
    // Read in parameters for the McGovern and De Pablo boundary consistent hills.
    } else if(!strcmp(word[0],"MCGDP_HILLS")){
      int read_lower_bound=0;
      int read_upper_bound=0;
      logical.mcgdp_hills=1;
      // first we select the proper CV
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv);}
      else{plumed_error("WITH MCGDP_HILLS YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");}
      hills.mcgdp_reshape_flag[icv - 1] = 1;
      // then we parse the line
      logical.upper[icv-1] = 1;
      logical.do_walls = 1;   // ### For modified output format
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")) {
          iw++;  // already read
        } else if(!strcmp(word[iw],"UPPER_BOUND")) {
          iw++; sscanf(word[iw], "%lf", &uno); hills.hill_upper_bounds[icv-1]=(real)uno; read_upper_bound=1;
        } else if(!strcmp(word[iw],"LOWER_BOUND")) {
          iw++; sscanf(word[iw], "%lf", &uno); hills.hill_lower_bounds[icv-1]=(real)uno; read_lower_bound=1;
        } else {
          plumed_error("Unknown flag for keyword MCGDP_HILLS");
        };
      }
      if(!read_lower_bound)
        plumed_error("WITH MCGDP_HILLS YOU ALWAYS HAVE TO SPECIFY THE \"LOWER_BOUND\" KEYWORD\n");
      if(!read_upper_bound)
        plumed_error("WITH MCGDP_HILLS YOU ALWAYS HAVE TO SPECIFY THE \"UPPER_BOUND\" KEYWORD\n");
      fprintf(mtd_data->fplog, "|-MCGDP_HILLS CONDITION %i: UPPER BOUND = %f, LOWER BOUND = %f \n\n",
             icv, hills.hill_upper_bounds[icv-1], hills.hill_lower_bounds[icv-1]); 
    // <JFD
    } else if(!strcmp(word[0],"INTERVAL")){
      int read_lower_limit=0;
      int read_upper_limit=0;
// first we select the proper CV
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv);}
      else{plumed_error("WITH INTERWAL YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");}
// then we parse the line
      logical.interval[icv-1]=1;  
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")) {
          iw++;  // already read
        } else if(!strcmp(word[iw],"LOWER_LIMIT")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvint.lower_limit[icv-1]=(real)uno; read_lower_limit=1;
        } else if(!strcmp(word[iw],"UPPER_LIMIT")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvint.upper_limit[icv-1]=(real)uno; read_upper_limit=1;
        } else {
          plumed_error("Unknown flag for keyword INTERVAL");
        };
      }
      if(!read_lower_limit)
        plumed_error("LOWER_LIMIT is missing");
      if(!read_upper_limit)
        plumed_error("UPPER_LIMIT is missing");
      fprintf(mtd_data->fplog, "|-INTERVAL ON COLVAR %i: LOWER_LIMIT = %f, UPPER_LIMIT = %f \n\n",
             icv, cvint.lower_limit[icv-1], cvint.upper_limit[icv-1]);  
    } else if(!strcmp(word[0],"STEER")){
      int read_max=0;
      int read_delta=0;
      int read_kappa=0;
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv); cvsteer.impose_start[icv-1] = 0;}
      else {plumed_error("WITH STEER YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");} 
      for(iw=1;iw<nw;iw++){ 
        if(!strcmp(word[iw],"CV")){
          iw++; // already read
        } else if(!strcmp(word[iw],"FROM")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.start[icv-1]=(real)uno; cvsteer.impose_start[icv-1]=1;
        } else if(!strcmp(word[iw],"TO")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.max[icv-1]=(real)uno; read_max=1;
        } else if(!strcmp(word[iw],"VEL")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.delta[icv-1]=(real) fabs(uno); read_delta=1;
        } else if(!strcmp(word[iw],"KAPPA")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.spring[icv-1]=(real) uno*mtd_data->eunit; read_kappa=1;
        } else if(!strcmp(word[iw],"RESTART")){
          logical.append = 1;
        } else {
          plumed_error("Unknown flag for keyword STEER");
        }
      }
      if(!read_max)  plumed_error("WITH STEER YOU ALWAYS HAVE TO SPECIFY THE \"TO\" KEYWORD\n");
      if(!read_delta)plumed_error("WITH STEER YOU ALWAYS HAVE TO SPECIFY THE \"VEL\" KEYWORD\n");
      if(!read_kappa)plumed_error("WITH STEER YOU ALWAYS HAVE TO SPECIFY THE \"KAPPA\" KEYWORD\n");

      logical.steer[icv-1]  = 1; 
      if(cvsteer.impose_start[icv-1]==0){
        fprintf(mtd_data->fplog, "|-STEERING COLVAR %i TO %f: VELOCITY=%lf cvunit/kstep, SPRING=%lf\n\n", icv,cvsteer.max[icv-1],cvsteer.delta[icv-1],cvsteer.spring[icv-1]/mtd_data->eunit);
      } else { 
        fprintf(mtd_data->fplog, "|-STEERING COLVAR %i FROM %f TO %f: VELOCITY=%lf cvunit/kstep, SPRING=%lf\n\n", icv, cvsteer.start[icv-1],cvsteer.max[icv-1],cvsteer.delta[icv-1],cvsteer.spring[icv-1]/mtd_data->eunit);
      } 
      if(logical.append) fprintf(mtd_data->fplog,"|-RESTARTING STEERING!\n");

  // DAFED  ########### ---------------------------------------------------------
  // The dafed structure is different from the steer structure
  // steer is a structure of arrays, dafed is an array of structures
    } else if(!strcmp(word[0],"DAFED")){
  #if defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
    	  int read_temperature=0;
		  int read_mass=0;
		  int read_kappa=0;
		  int read_tauthermo=0;
		  int read_nrespa_ggmt=0;
		  int n_respa_ggmt;

		  iw = seek_word(word,"CV");
		  if(iw>=0){ sscanf(word[iw+1], "%i", &icv);}
		  else {plumed_error("WITH DAFED YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");}

                dafed[icv-1].do_jacobian_force = 0;
                dafed[icv-1].do_periodicity = 0;
                dafed[icv-1].periodicity_low = 0.0;
                dafed[icv-1].periodicity_high = 0.0;
                dafed[icv-1].periodicity_gap = 0.0;

		  for(iw=1;iw<nw;iw++){
				if(!strcmp(word[iw],"CV")){
				  iw++; // already read
				} else if(!strcmp(word[iw],"TEMPERATURE")) {
				  iw++; sscanf(word[iw], "%lf", &uno); dafed[icv-1].temperature=(real)uno; read_temperature=1;
				} else if(!strcmp(word[iw],"MASS")) {
				  iw++; sscanf(word[iw], "%lf", &uno); dafed[icv-1].mass=(real)uno; read_mass=1;
				} else if(!strcmp(word[iw],"KAPPA")) {
				  iw++; sscanf(word[iw], "%lf", &uno); dafed[icv-1].kappa=(real)uno*mtd_data->eunit; read_kappa=1;
				} else if(!strcmp(word[iw],"TAUTHERMO")) {
				  iw++; sscanf(word[iw], "%lf", &uno); dafed[icv-1].tauthermo=(real)uno; read_tauthermo=1;
				} else if(!strcmp(word[iw],"N_RESPA_GGMT")) {
				  iw++; sscanf(word[iw], "%lf", &uno); dafed[icv-1].ggmt.n_respa_ggmt=(int)uno; read_nrespa_ggmt=1;
                                } else if(!strcmp(word[iw],"PERIODIC")) {
                                  iw++;
				  if (!strcmp(word[iw],"MINUS_PI")) {
					dafed[icv-1].periodicity_low= - M_PI;
				  }else{
					sscanf(word[iw], "%lf", &uno);
					dafed[icv-1].periodicity_low=(real)uno;
				  }
				  iw++;
				  if (!strcmp(word[iw],"PLUS_PI")) {
				  dafed[icv-1].periodicity_high= + M_PI;
				  }else if (!strcmp(word[iw],"PLUS_2PI")){
				  dafed[icv-1].periodicity_high= + 2.0* M_PI;
				  }else{
				  sscanf(word[iw], "%lf", &uno);
				  dafed[icv-1].periodicity_high=(real)uno;
				}
				dafed[icv-1].periodicity_gap = dafed[icv-1].periodicity_high - dafed[icv-1].periodicity_low;
				dafed[icv-1].do_periodicity = 1;
				} else if(!strcmp(word[iw],"JACOBIAN_FORCE")) {
				  dafed[icv-1].do_jacobian_force = 1;
				} else {
				  plumed_error("Unknown flag for keyword DAFED");
				}
		}
		if(!read_temperature)  plumed_error("WITH DAFED YOU ALWAYS HAVE TO SPECIFY THE \"TEMPERATURE\" KEYWORD\n");
		if(!read_mass)  plumed_error("WITH DAFED YOU ALWAYS HAVE TO SPECIFY THE \"MASS\" KEYWORD\n");
		if(!read_kappa)  plumed_error("WITH DAFED YOU ALWAYS HAVE TO SPECIFY THE \"KAPPA\" KEYWORD\n");
		if(!read_tauthermo)  plumed_error("WITH DAFED YOU ALWAYS HAVE TO SPECIFY THE \"TAUTHERMO\" KEYWORD\n");

		// default values
		if(!read_nrespa_ggmt) dafed[icv-1].ggmt.n_respa_ggmt=1;

		logical.dafed[icv-1]  = 1;
		logical.do_dafed  = 1;

		fprintf(mtd_data->fplog, "|- DAFED ON COLVAR %i WITH FOLLOWING PARAMETERS\n",icv);
		fprintf(mtd_data->fplog, "|-\tTEMPERATURE %f\n",dafed[icv-1].temperature);
		fprintf(mtd_data->fplog, "|-\tMASS %f\n",dafed[icv-1].mass);
		fprintf(mtd_data->fplog, "|-\tKAPPA %f\n",dafed[icv-1].kappa);
		fprintf(mtd_data->fplog, "|-\tTHERMOSTAT TAU %f\n",dafed[icv-1].tauthermo);
		if (dafed[icv-1].do_jacobian_force) {
			fprintf(mtd_data->fplog, "|-\tApply JACOBIAN_FORCE\n");
		}

  #else
	plumed_error("DAFED is not implemented with MD packages other than Gromacs 4.x");
  #endif
    } else if(!strcmp(word[0],"DAFED_CONTROL")){

    	int read_nrespa=0;
    	for(iw=1;iw<nw;iw++){
    			if(!strcmp(word[iw],"RESTART")){
    					  dafed_control.restart = 1;
    					  iw++; sscanf(word[iw], "%s", &(dafed_control.in_file));
    					  fprintf(mtd_data->fplog, "|- WILL RESTART WITH d-AFED STATE FROM FILE %s \n",dafed_control.in_file );
    			} else if(!strcmp(word[iw],"WRITE_STATE")) {
    					  iw++; sscanf(word[iw], "%lf", &uno); dafed_control.write_freq= (int)uno;
    					  fprintf(mtd_data->fplog, "|- WILL WRITE d-AFED STATE TO FILE DAFED_STATE EVERY %d STEPS \n",dafed_control.write_freq);
    			} else if(!strcmp(word[iw],"N_RESPA")) {
						  iw++; sscanf(word[iw], "%lf", &uno); dafed_control.n_respa=(int)uno; read_nrespa=1;
						  fprintf(mtd_data->fplog, "|- RESPA STEPS FOR d-AFED %d\n",dafed_control.n_respa);
    			}
    	}
    	strcpy(dafed_control.out_file,"DAFED_STATE");  // Default output file name
		if(!read_nrespa) dafed_control.n_respa=1;

  // #### -------------------------------------------------------------------------------

    } else if(!strcmp(word[0],"UMBRELLA")){
      logical.do_walls = 1;   // ### For modified output format
      int read_kappa=0;
      int read_slope=0;
      int read_anneal=0;
      int read_at=0;
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv); cvsteer.impose_start[icv-1]=1; cvsteer.delta[icv-1]=0.;}
      else {plumed_error("WITH UMBRELLA YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");} 
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")){
          iw++; // already read
        } else if(!strcmp(word[iw],"AT")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.max[icv-1]=cvsteer.start[icv-1]=(real) uno; read_at=1;
        } else if(!strcmp(word[iw],"KAPPA")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.spring[icv-1]=(real) uno*mtd_data->eunit; read_kappa=1;
        } else if(!strcmp(word[iw],"SLOPE")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.slope[icv-1]=(real) uno*mtd_data->eunit; read_slope=1;
#if defined (PLUMED_GROMACS4)|| defined (PLUMED_GROMACS45)
        } else if(!strcmp(word[iw],"ANNEALING")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.annealing[icv-1]=(real) uno*mtd_data->eunit; read_anneal=1;
#endif
        } else if(!strcmp(word[iw],"RESTART")){
          logical.append = 1;
        } else {
          plumed_error("Unknown flag for keyword UMBRELLA");
        }
      }
      if(!(read_kappa||read_slope))plumed_error("WITH UMBRELLA YOU ALWAYS HAVE TO SPECIFY THE \"KAPPA\" OR THE \"SLOPE\" KEYWORD\n");
      if(!read_at)plumed_error("WITH UMBRELLA YOU ALWAYS HAVE TO SPECIFY THE \"AT\" KEYWORD\n");

      logical.steer[icv-1]  = 1; 
      fprintf(mtd_data->fplog, "|-UMBRELLA SAMPLING OF COLVAR %i AT %f: SPRING=%lf SLOPE=%lf\n\n",
              icv,cvsteer.max[icv-1],cvsteer.spring[icv-1]/mtd_data->eunit,cvsteer.slope[icv-1]/mtd_data->eunit);
      if(logical.append) fprintf(mtd_data->fplog,"|-RESTARTING UMBRELLA SAMPLING!\n");
      if(read_anneal) fprintf(mtd_data->fplog,"|-UMBRELLA %i DEPEND ON TEMPERATURE AS %lf/T\n",icv-1, cvsteer.annealing[icv-1]/mtd_data->eunit);

    } else if(!strcmp(word[0],"CONSTRAINT")){
      int read_at=0;
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv); cvcnstr.delta[icv-1]=1.e-5 ; cvcnstr.maxiter[icv-1]=1000 ;  cvcnstr.spring[icv-1]=1.; cvcnstr.verbose[icv-1]=0;}
      else {plumed_error("WITH CONSTRAINT YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");} 
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")){
          iw++; // already read
        } else if(!strcmp(word[iw],"AT")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvcnstr.pos[icv-1]=(real) uno; read_at=1;
        } else if(!strcmp(word[iw],"KAPPA")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvcnstr.spring[icv-1]=(real) uno*mtd_data->eunit; 
        } else if(!strcmp(word[iw],"DELTA")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvcnstr.delta[icv-1]=(real) uno; 
        } else if(!strcmp(word[iw],"MAXITER")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvcnstr.maxiter[icv-1]=(real) uno; 
        } else if(!strcmp(word[iw],"VERBOSE")) {
          cvcnstr.verbose[icv-1]=1; 
        } else {
          plumed_error("Unknown flag for keyword CONSTRAINT");
        }
      }
      if(!read_at)plumed_error("WITH CONSTRAINT YOU ALWAYS HAVE TO SPECIFY THE \"AT\" KEYWORD\n");

      logical.cnstr[icv-1]  = 1; 
      fprintf(mtd_data->fplog, "|-CONSTRAINED SAMPLING OF COLVAR %i AT %f: SPRING=%lf\n\n", icv,cvcnstr.pos[icv-1],cvcnstr.spring[icv-1]/mtd_data->eunit);
      logical.do_constraint=1; 
    } else if(!strcmp(word[0],"STEERPLAN")){
         if(logical.do_steerplan==1)
              plumed_error("ONLY ONE STEERPLAN IS ALLOWED...BUT YOU MAY HAVE MANY CVS IN THERE!!!");
         logical.do_steerplan=1;    
         read_steerplan(word, count, &input,  &iline ,mtd_data->fplog);
    } else if(!strcmp(word[0],"ABMD")){
      int read_max=0;
      int read_from=0;
      int read_kappa=0;
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv); }
      else {plumed_error("WITH ABMD YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");} 
      for(iw=1;iw<nw;iw++){ 
        if(!strcmp(word[iw],"CV")){
          iw++; // already read
        } else if(!strcmp(word[iw],"RESTART")) {
          iw++; sscanf(word[iw], "%lf", &uno); abmd.min[icv-1]=(real)uno; read_from=1; logical.restart_abmd=1;
        } else if(!strcmp(word[iw],"TO")) {
          iw++; sscanf(word[iw], "%lf", &uno); abmd.exp[icv-1]=(real)uno; read_max=1;
        } else if(!strcmp(word[iw],"KAPPA")) {
          iw++; sscanf(word[iw], "%lf", &uno); abmd.spring[icv-1]=(real) uno*mtd_data->eunit; read_kappa=1;
        } else {
          plumed_error("Unknown flag for keyword ABMD");
        }
      }
      if(!read_max)  plumed_error("WITH ABMD YOU MUST SPECIFY THE \"TO\" KEYWORD\n");
      if(!read_kappa)plumed_error("WITH ABMD YOU MUST SPECIFY THE \"KAPPA\" KEYWORD\n");

      logical.abmd[icv-1]  = 1; 
      if(read_from==0){
        abmd.min[icv-1] = 9999999999.0;
        fprintf(mtd_data->fplog, "|-ABMD ON COLVAR %i TO %f: , SPRING=%lf\n\n", icv, abmd.exp[icv-1], abmd.spring[icv-1]/mtd_data->eunit);
      } else { 
        fprintf(mtd_data->fplog, "|-ABMD ON COLVAR %i RESTARTING FROM %f TO %f: , SPRING=%lf\n\n", icv, abmd.min[icv-1], abmd.exp[icv-1], abmd.spring[icv-1]/mtd_data->eunit);
      } 
      if(logical.restart_abmd) fprintf(mtd_data->fplog,"|-RESTARTING ABMD!\n");

    } else if(!strcmp(word[0],"NOHILLS")){
      iw = seek_word(word,"CV");
      if(iw>=0) sscanf(word[iw+1], "%i", &icv);
      colvar.on[icv-1] = 0;
      fprintf(mtd_data->fplog, "|-NO HILLS ON COLVAR %i\n", icv);
    } else if(!strcmp(word[0],"INVERT")){
      logical.do_inversion = 1;
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv);
      } else{plumed_error("WITH INVERT YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");}
      fprintf(mtd_data->fplog, "APPLY INVERSION TO FREE ENERGY ON CV %i:\n", icv);
      iw=seek_word(word,"REFLECTION");
      if (iw>=0) { sscanf(word[iw+1],"%lf",&uno);
        colvar.inv_ref[icv-1]= (real) uno;
      }
      iw=seek_word(word,"INVERSION");
      if (iw>=0) { sscanf(word[iw+1],"%lf",&uno);
        colvar.inv_inv[icv-1]= (real) uno;
      }
      if (colvar.inv_inv[icv-1]<colvar.inv_ref[icv-1]) plumed_error("REFLECTION INTERVAL LARGER THAN THE INVERSION ONE");       
      iw=seek_word(word,"MAXHEIGHT");
      if (iw>=0) { sscanf(word[iw+1],"%lf",&uno);
        colvar.inv_maxww[icv-1]= (real) uno;
      }
      iw=seek_word(word,"LIMIT1");
      if (iw>=0) { sscanf(word[iw+1],"%lf",&uno);
        colvar.inv_limit[icv-1][0]= (real) uno;
        logical.invert[icv-1][0] = 1;
        fprintf(mtd_data->fplog, " LIMIT1 = %f\n",colvar.inv_limit[icv-1][0]); }
      iw=seek_word(word,"LIMIT2");
      if (iw>=0) { sscanf(word[iw+1],"%lf",&uno);
        colvar.inv_limit[icv-1][1]= (real) uno;
        logical.invert[icv-1][1] = 1;
        fprintf(mtd_data->fplog, " LIMIT2 = %f\n",colvar.inv_limit[icv-1][1]); }
      if (logical.invert[icv-1][0]==0) {
        if (logical.invert[icv-1][1]==0) plumed_error("NO LIMITS FOR INVERSION FOUND: SPECIFY AT LIST ONE!");
      }
      fprintf(mtd_data->fplog, " REFLECTION INTERVAL (gaussian width units) = %f\n",colvar.inv_ref[icv-1]);
      fprintf(mtd_data->fplog, " INVERSION INTERVAL (gaussian width units) = %f\n",colvar.inv_inv[icv-1]);
      fprintf(mtd_data->fplog, " MAX GAUSSIAN HEIGHT FACTOR = %f\n\n",colvar.inv_maxww[icv-1]); 
    } else if(!strcmp(word[0],"UREFLECT")){
      iw = seek_word(word,"CV");
      if(iw>=0) sscanf(word[iw+1], "%i", &icv);
      iw = seek_word(word,"LIMIT");
      if(iw>=0) sscanf(word[iw+1], "%lf", &uno);
      cvw.upper[icv-1] = (real) uno;
      logical.ureflect[icv-1] = 1;
      fprintf(mtd_data->fplog, "|-UPPER REFLECTING WALL ON CV %i, AT %f\n\n", icv, cvw.upper[icv-1]);
    } else if(!strcmp(word[0],"LREFLECT")){
      iw = seek_word(word,"CV");
      if(iw>=0) sscanf(word[iw+1], "%i", &icv);
      iw = seek_word(word,"LIMIT");
      if(iw>=0) sscanf(word[iw+1], "%lf", &uno);
      cvw.lower[icv-1] = (real) uno;
      logical.lreflect[icv-1] = 1;
      fprintf(mtd_data->fplog, "|-LOWER REFLECTING WALL ON CV %i, AT %f\n\n", icv, cvw.lower[icv-1]);
    } else if(!strcmp(word[0],"DEBUG_GRID")){
      logical.debug_grid=1;
    } else if(!strcmp(word[0],"NOSPLINE")){
      logical.donot_spline=1;
      fprintf(mtd_data->fplog, "|- GRID SPLINE TURNED OFF\n");
    } else if(!strcmp(word[0],"GRID")){
      int read_cv=0;
      int read_min=0;
      int read_max=0;
      int read_nbin=0;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")){
          iw++; sscanf(word[iw], "%i", &icv);  read_cv=1;
        } else if(!strcmp(word[iw],"MIN")){
          iw++; sscanf(word[iw], "%s", tmpmeta); uno = plumed_atof(tmpmeta); bias_grid.min[bias_grid.ncv] = (real) uno; read_min=1;
        } else if(!strcmp(word[iw],"MAX")){
          iw++; sscanf(word[iw], "%s", tmpmeta); due = plumed_atof(tmpmeta); bias_grid.max[bias_grid.ncv] = (real) due; read_max=1;
        } else if(!strcmp(word[iw],"NBIN")){
          iw++; sscanf(word[iw], "%i", &(bias_grid.bin[bias_grid.ncv])); read_nbin=1;
        } else if(!strcmp(word[iw],"PBC")){
          bias_grid.period[bias_grid.ncv] = 1;
        } else {
          plumed_error("Unknown flag for keyword GRID");
        }
      };
      if(!read_cv) plumed_error("WITH GRID YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");
      if(!read_min) plumed_error("WITH GRID YOU ALWAYS HAVE TO SPECIFY THE \"MIN\" KEYWORD\n");
      if(!read_max) plumed_error("WITH GRID YOU ALWAYS HAVE TO SPECIFY THE \"MAX\" KEYWORD\n");
      if(!read_nbin) plumed_error("WITH GRID YOU ALWAYS HAVE TO SPECIFY THE \"NBIN\" KEYWORD\n");
      fprintf(mtd_data->fplog, "|-GRID ACTIVE ON CV %i NBIN %d MIN %f MAX %f \n", icv, bias_grid.bin[bias_grid.ncv],bias_grid.min[bias_grid.ncv],bias_grid.max[bias_grid.ncv]);
      if(bias_grid.period[bias_grid.ncv]) fprintf(mtd_data->fplog, "|-- PERIODIC GRID IS ON\n");
      bias_grid.index[bias_grid.ncv] = icv-1;
      for(i=0;i<bias_grid.ncv;i++) if(bias_grid.index[i]==bias_grid.index[bias_grid.ncv]) plumed_error("GRID is already ACTIVE for this CV");
      bias_grid.ncv += 1;
      logical.do_grid = 1;
      fprintf(mtd_data->fplog, "\n");
    } else if(!strcmp(word[0],"WRITE_GRID")){
      logical.write_grid=1;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"W_STRIDE")){
          iw++; sscanf(word[iw], "%i", &(bias_grid.w_stride));
        } else if(!strcmp(word[iw],"FILENAME")){
          iw++; sscanf(word[iw], "%s", bias_grid.w_file);
        }
      };
      fprintf(mtd_data->fplog, "|-WRITING GRID ON FILE \n");
      fprintf(mtd_data->fplog, "|--STRIDE %d FILENAME %s\n\n", bias_grid.w_stride, bias_grid.w_file);
    } else if(!strcmp(word[0],"READ_GRID")){
      logical.read_grid=1;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"FILENAME")){
          iw++; sscanf(word[iw], "%s", bias_grid.r_file);
        }
      };  
      fprintf(mtd_data->fplog, "|-READING GRID FROM FILE %s \n\n",bias_grid.r_file);
    } else if(!strcmp(word[0],"PROJ_GRAD")){ 
       iw = seek_word(word,"CV"); // look for a group
       if (iw>=0){
          colvar.pg.nlist=plumed_get_group(word[iw+1],&colvar.pg.list,0,&input,mtd_data->fplog); 
       } else {
           plumed_error("WITH PROJ_GRAD YOU ALWAYS HAVE  TO SPECIFY THE \"CV\" KEYWORD\n AND ASSOCIATE A GROUP TO THIS \n");
       }
  //     iw = seek_word(word,"DUMPFILE");
  //     if (iw>=0){
  //        sscanf(word[iw+1], "%s", pg.dumpfile);
  //     } else {
  //        strcpy(pg.dumpfile,"proj_grad.dat")
  //     }
  //     iw = seek_word(word,"W_STRIDE"); 
  //     if (iw>=0){
  //        sscanf(word[iw+1], "%lf", &uno);pg.w_stride=(real)uno;
  //     } else {
  //       pg.w_stride=1.0;
  //     }
    } else if(!strcmp(word[0],"EXTERNAL")){
      if(logical.do_external) plumed_error("ONLY ONE EXTERNAL POTENTIAL ALLOWED");
      logical.do_external=1;
      for(iw=1;iw<nw;iw++) if(!strcmp(word[iw],"NCV")){iw++; sscanf(word[iw], "%d", &(extpot.ncv));}
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")){
          iw++; for(i=0;i<extpot.ncv;i++) {sscanf(word[iw], "%d", &(extpot.index[i])); extpot.index[i] -= 1; iw++;}; 
          iw--;
        } else if(!strcmp(word[iw],"FILENAME")){
          iw++; sscanf(word[iw], "%s", extpot.r_file);
        }
      };
      fprintf(mtd_data->fplog, "|-EXTERNAL POTENTIAL ON CV ");
      for(i=0;i<extpot.ncv;i++) fprintf(mtd_data->fplog, " %d ",extpot.index[i]+1);
      fprintf(mtd_data->fplog, "FROM FILE %s \n\n",extpot.r_file);
    } else if(!strcmp(word[0],"STOPWHEN")) {
#if ! defined(PLUMED_GROMACS45)
      fprintf(mtd_data->fplog, "WARNING: You are not using GROMACS 4.5\n");
      fprintf(mtd_data->fplog, "WARNING: STOPWHEN will not dump the final configuration properly\n");
#endif
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv);icv--;
         fprintf(mtd_data->fplog, "|-STOPWHEN ENABLED ON CV %d: \n",icv+1);
      }
      else{plumed_error("WITH STOPWHEN YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");}  
      iw = seek_word(word,"MORETHAN");
      if(iw>=0){ sscanf(word[iw+1], "%lf",&uno);stopwhen.max[icv]=(real)uno; stopwhen.actmax[icv]=1; 
         fprintf(mtd_data->fplog, "|-STOPWHEN CV %d IS MORE THAN %lf\n",icv+1,stopwhen.max[icv]);
      }
      iw = seek_word(word,"LESSTHAN");
      if(iw>=0){ sscanf(word[iw+1], "%lf",&uno);stopwhen.min[icv]=(real)uno; stopwhen.actmin[icv]=1; 
         fprintf(mtd_data->fplog, "|-STOPWHEN CV %d IS LESS THAN %lf\n",icv+1,stopwhen.min[icv]);
      }
      if (stopwhen.actmin[icv]==0 && stopwhen.actmax[icv] ==0){
           fprintf(mtd_data->fplog, "|-STOPWHEN SYNTAX:\n" );
           fprintf(mtd_data->fplog, "|-    STOPWHEN CV 1  MORETHAN 3.0 LESSTHAN 1.0 \n" );
           fprintf(mtd_data->fplog, "|-    STOPWHEN CV 1  LESSTHAN 1.0 \n" );
           fprintf(mtd_data->fplog, "|-     \n" );
           char buf[1024];
           sprintf(buf, "STOPWHEN DIED WITH ERRORS: CHECK THE INPUT!!!!!\n");
           plumed_error(buf);
      };
      fprintf(mtd_data->fplog, "\n");
    } else if(!strcmp(word[0],"COUPLINGMATRIX")) { 
		read_couplingmatrix  ( word, &input, mtd_data->fplog );
	  } else {
      char buf[1024];
      sprintf(buf, "Line %i Unkwown Keyword %s \n", iline+1, word[0]);
      plumed_error(buf);
    }
  }
// clean input parser
  plumed_clear_input(&input);

// set the number of collective variables
  colvar.nconst = count;

// Call reconnaissance metadynamics setup
#ifdef RECONMETAD 
#ifndef DRIVER
   if( reconOn==1 && colvar.nbespoke>0 ){
      plumed_error("can't do reconnaissance metadynamics with bespoke collective coordinates");
   } else if( reconOn==1 ){ 
      // setup a run in which all cvs are used in reconnaissance if not otherwise instructed in input
      if( reconinpt.nconst==0 ){ 
        reconinpt.nconst=colvar.nconst; srenew(reconinpt.cvlist, reconinpt.nconst); 
        for(i=0;i<reconinpt.nconst;++i){ reconinpt.cvlist[i]=i; } 
      }   

      // transfer the periods 
      double periods[reconinpt.nconst]; 
      for(i=0;i<reconinpt.nconst;i++){
         periods[i]=0.;
         if( colvar.type_s[reconinpt.cvlist[i]]==5 ){
            if( colvar.doTrig[reconinpt.cvlist[i]]==0 ){ periods[i]=2.*M_PI; }
            else{ periods[reconinpt.cvlist[i]]=0.; }
         }   
      }    
 
      // and create the reconnaissance metadynamics object
      create_recon(&myreconObj);
      double tstep; tstep=mtd_data->dt;
      setup_recon( periods, tstep, reconinpt, myreconObj, mtd_data->fplog);

      fprintf(mtd_data->fplog, "|- RECONNAISSANCE METADYNAMICS ON COLVARS :");
      for(i=0;i<reconinpt.nconst;i++){
        fprintf(mtd_data->fplog," %d ",reconinpt.cvlist[i]+1); if((i+1)%20==0)fprintf(mtd_data->fplog,"\n                    ");
      }
      fprintf(mtd_data->fplog,"\n\n");
   
   // This is the setup for bespoke collective coordinates ( this might have to be changed in the future so that it works more like the above )
   } else if( colvar.nbespoke>0 ){
#else
   if( colvar.nbespoke>0 ){
#endif
      int ncolvar; ncolvar=colvar.nconst-colvar.nbespoke;
      for(i=ncolvar;i<colvar.nconst;i++){
        if( colvar.type_s[i]!=46){ plumed_error("BESPOKE COLLECTIVE COORDINATES MUST COME AFTER ALL OTHER COLLECTIVE COORDINATES"); }
      }
      double periods[ncolvar];
      for(i=0;i<ncolvar;i++){
         periods[i]=0.;
         if( colvar.type_s[i]==5 ){
            if( colvar.doTrig[i]==0 ){ periods[i]=2.*M_PI; }
            else{ periods[i]=0.; }
         }   
         if( colvar.type_s[i]==46){ plumed_error("BESPOKE COLLECTIVE COORDINATES MUST COME AFTER ALL OTHER COLLECTIVE COORDINATES"); }
      }    

#ifdef CVS
      create_bespoke( ncolvar, colvar.nbespoke, &mybespokeObj);
      setup_bespoke( bespoke_input, ncolvar, colvar.nbespoke, periods, mybespokeObj, mtd_data->fplog );
#endif
   }
#endif

// unset the SIGMA<0 CVs for hills
  for(i=0;i<colvar.nconst;i++){
     if(colvar.delta_r[i]<0.){colvar.on[i]=0;}
     if(logical.do_hills){
       if(colvar.on[i]){fprintf(mtd_data->fplog, "|-HILLS ACTIVE ON COLVAR %i\n", i+1);} 
       else            {fprintf(mtd_data->fplog, "|-NO HILLS     ON COLVAR %i\n", i+1);}
     }else if(logical.tamd){
       if(colvar.on[i]){
         tamd.spring[i]=mtd_data->boltz*tamd.simtemp/(colvar.delta_r[i]*colvar.delta_r[i]);
         fprintf(mtd_data->fplog, "|-CV %i : TAMD/DAFED WITH SPRING CONSTANT = %lf\n",i+1,tamd.spring[i]);
       } 
       else{
         fprintf(mtd_data->fplog, "|-CV %i : NO TAMD/DAFED\n", i+1);
       }
     }
  }

// check the correctenes of the input parsed
  if(colvar.nconst > nconst_max) {
    plumed_error("Too many colvars. Change NCONST_MAX in metadyn.h !!!!!!!!!!!\n");	
  }

// checking for conflicts in directive keywords
#ifdef RECONMETAD  
  if(!logical.do_hills && !logical.commit && !logical.do_dafed && reconOn!=1){
    //if( reconOn==1 && reconinpt.monitor==1 ){
    //   fprintf(mtd_data->fplog,"|-ANALYSIS: YOU WILL MONITOR BASIN OCCUPANCIES ONLY\n\n");  
    //} 
    // if( reconOn!=1 ){
    fprintf(mtd_data->fplog, "|-ANALYSIS: YOU WILL ONLY MONITOR YOUR CVs DYNAMICS\n\n");
    // }
  }
#else
  if(!logical.do_hills && !logical.commit && !logical.do_dafed ){   // ### DAFED
    fprintf(mtd_data->fplog, "|-ANALYSIS: YOU WILL ONLY MONITOR YOUR CVs DYNAMICS\n\n");
  }
#endif

// derivatives debug with ENERGY CV not allowed
  if(logical.energy && (logical.debug || logical.debug_derivatives)) 
    plumed_error("DERIVATIVES DEBUG with ENERGY CV not allowed");

  if(logical.welltemp && !logical.do_hills)  plumed_error("WELLTEMPERED must be used with HILLS keyword");
  // JFD>
  if(logical.dicksonian_tempering && !logical.do_grid)  plumed_error("DICKSONIAN_TEMPERING must be used with GRID keyword");
  if(logical.transition_tempering && !logical.do_hills) plumed_error("TRANSITIONTEMPERED must be used with HILLS keyword");
  if(logical.transition_tempering && !logical.do_grid)  plumed_error("TRANSITIONTEMPERED must be used with GRID keyword");
  if(logical.ttdebug && !logical.transition_tempering)  plumed_error("DEBUG_TRANSITIONTEMPERED must be used with TRANSITIONTEMPERED keyword");
  if(logical.mcgdp_hills && !logical.do_hills)  plumed_error("MCGDP_HILLS must be used with HILLS keyword");
  if(logical.mcgdp_hills && !logical.do_grid)  plumed_error("MCGDP_HILLS must be used with GRID keyword");
  // <JFD

  if(logical.commit && logical.do_hills) plumed_error("KEYWORD 'COMMITMENT' AND 'HILLS' ARE NOT COMPATIBLE");

// in case of parallel rescale hills heigth with temperature
  if(logical.do_hills&&logical.remd&&(!logical.rpxm)&&(!logical.norescale)) hills.wwr *= mtd_data->rteio/mtd_data->rte0;

// in case of PTMETAD and well-tempered set the right simtemp
  if(logical.do_hills&&logical.remd&&(!logical.rpxm)&&logical.welltemp) colvar.simtemp = mtd_data->rteio;

// check for untested features
  if(!logical.enable_untested_features) {
   if(logical.debug_derivatives) plumed_error("DEBUG_DERIVATIVES NOT ENABLED");
   if(colvar.ptmetad_neighbours) plumed_error("NEIGHBOUR HILLS NOT ENABLED");
   if(logical.debug_grid)        plumed_error("DEBUG_GRID NOT ENABLED");
   if(hills.max_height>0.0)      plumed_error("MAX_HEIGHT NOT ENABLED");
   if(logical.hrex)              plumed_error("HAMILTONIAN REPLICA-EXCHANGE NOT ENABLED");
   if(logical.tamd)              plumed_error("TAMD/DAFED NOT ENABLED");
   // JFD>
   if(logical.dicksonian_tempering)  plumed_error("WELLTEMPERED USING THE DICKSON RULE NOT ENABLED");
   if(logical.transition_tempering)  plumed_error("TRANSITIONTEMPERED NOT ENABLED");
   if(logical.mcgdp_hills)  plumed_error("MCGDP_HILLS NOT ENABLED");
   // <JFD
  }

// GRID and WRITE/READ
  if(logical.read_grid  && !logical.do_grid) plumed_error("GRID must be active to use READ_GRID\n");
  if(logical.write_grid && !logical.do_grid) plumed_error("GRID must be active to use WRITE_GRID\n");
  if(logical.do_walkers &&  logical.read_grid) plumed_error("READ_GRID cannot be used with MULTIPLE_WALKERS\n");

// checking if GRID and HILLS active variables are consistent
  if(logical.do_grid) {
   icv = 0;
   for(i=0;i<colvar.nconst;i++) if(colvar.on[i]) icv++;  
   if(icv!=bias_grid.ncv) plumed_error("Inconsistency between GRID and HILLS variables. Please, check !!!!!!!!!!!\n"); 
   for(i=0;i<bias_grid.ncv;i++) if(!colvar.on[bias_grid.index[i]] || bias_grid.index[i]>=colvar.nconst) 
     plumed_error("Inconsistency between GRID and HILLS variables. Please, check !!!!!!!!!!!\n");
// in case initialize grid stuff
   grid_initialize(&bias_grid);
  }
 
// check EXTERNAL potential CVs
 if(logical.do_external){
  if(extpot.ncv>colvar.nconst) plumed_error("Too many CVs for EXTERNAL potential.\n");
  for(i=0;i<extpot.ncv;i++) if(extpot.index[i]>=colvar.nconst) plumed_error("Check the CVs for EXTERNAL potential. Do they exist?\n");
 }

// multiple walkers allocation
  hills.line_counter = (fpos_t *)calloc(hills.nwalkers,sizeof(fpos_t));

// check for needed projection
  if(colvar.pg.nlist!=0){
        // make the projection tables
       fprintf(mtd_data->fplog, "|- FOUND PROJ_GRAD KEYWORD: NCV involved %d\n",colvar.pg.nlist);
       int j; 
       fprintf(mtd_data->fplog, "|- WHICH ARE: ");
       for(j=0;j<colvar.pg.nlist;j++){fprintf(mtd_data->fplog, " %d",colvar.pg.list[j]);}
       fprintf(mtd_data->fplog, "\n");
       setup_projections( &(colvar.pg));          
  }
// check whether a variable has to be calculated every step or not
  fprintf(mtd_data->fplog, "|- DIFFERENT COLLECTIVE VARIABLE WILL BE CALCULATED AT DIFFERENT TIMES\n");
  // set always parameter  
  // note colvar.on is for activating metadynamics
  // logical.always is for only calculating cvs  
  for(i=0;i<colvar.nconst;i++) {
    if(colvar.on[i]||logical.steer[i]||logical.abmd[i]||logical.upper[i]||logical.lower[i]||logical.cnstr[i]|| steerplan.isactive[i] ||logical.debug_derivatives || stopwhen.actmax[i] || stopwhen.actmin[i]) logical.always[i]=1;
  }

// check if INTERVAL is used more than once
  tmpc = 0;
  for(i=0;i<colvar.nconst;i++) if(logical.interval[i]) tmpc++;
  if(tmpc>1) plumed_error("INTERVAL CAN BE USED ONLY ON A SINGLE DIMENSION, ON MANY DIMENSIONS ITS BEHAVOUR IS NOT TESTED!\n");

// add histogram variables
  for(i=0;i<colvar.nconst;i++) { 
     if(colvar.type_s[i]==49){
        for(j=0;j<colvar.histo_ncv[i];j++){ logical.always[colvar.histo_cvlist[i][j]]=1; }
     } 
  }   
// add committment variables
  if(logical.commit)      for(i=0;i<commit.ncv;i++) logical.always[commit.index[i]]=1;
// add external potential
  if(logical.do_external) for(i=0;i<extpot.ncv;i++) logical.always[extpot.index[i]]=1;  
#ifdef RECONMETAD
  // and reconnaissance metadynamics
  if( reconOn==1 ){                   // && reconinpt.monitor!=1 ){
    for(i=0;i<reconinpt.nconst;++i){ logical.always[reconinpt.cvlist[i]]=1; }
  }
#endif
#ifdef CVS
  // and cvs from which bespoke cvs are calculated
  for(i=0;i<colvar.bespoke_ncv;i++){ logical.always[colvar.bespoke_cvlist[i]]=1; }
#endif

// #### d-AFED initialization ----------------------------------------
  for(i=0;i<colvar.nconst;i++){
		if (logical.dafed[i]) {
			initialize_dafed(&dafed[i], mtd_data->dt);
			if (dafed_control.restart) {
				dafed[i].do_initialize_s = 0;
				dafed[i].do_skip_integration = 1;
			}
			logical.always[i]=1;
		}
  }
  // This overwrites the default values just initialized above.
  if (dafed_control.restart){
 	  read_dafed_state();
  }
// ####----------------------------------------------------------------

  for(i=0;i<colvar.nconst;i++) {
    if(logical.always[i]) fprintf(mtd_data->fplog, "|--CV %2i WILL BE EVALUATED EACH STEP\n", i+1);
    else fprintf(mtd_data->fplog, "|--CV %2i WILL BE EVALUATED ONLY WHEN NEEDED (OUTPUT OR EXCHANGE TRIAL)\n", i+1);
  }


// printout PLEASE_CITE
  cite_please("bono+09cpc",mtd_data->fplog);
  if(logical.do_hills) cite_please("laio-parr02pnas",mtd_data->fplog);
  if(logical.remd && !logical.rpxm) cite_please("buss+06jacs",mtd_data->fplog);
  if(logical.rpxm)       cite_please("pian-laio07jpcb",mtd_data->fplog);
  if(logical.welltemp)   cite_please("bard+08prl",mtd_data->fplog);
  // JFD>
  if(logical.dicksonian_tempering) cite_please("dickson2011pre",mtd_data->fplog);
  if(logical.mcgdp_hills) cite_please("mcgovern2013jcp", mtd_data->fplog);
  // <JFD
  if(logical.path)       cite_please("bran+07jcp",mtd_data->fplog);
  if(logical.puckering)  cite_please("sega+09jcp",mtd_data->fplog);
  if(logical.do_walkers) cite_please("rait+06jpcb",mtd_data->fplog);
  if(logical.do_alphabetarmsd) cite_please("pietrucci+09jctc",mtd_data->fplog);
  if(logical.do_sprint) cite_please("pietrucci+11prl",mtd_data->fplog);
  if(logical.do_inversion) cite_please("marinell-crespo10",mtd_data->fplog);
  if (logical.do_dafed) { // d-AFED ####
	  cite_please("abrams08jpcb",mtd_data->fplog);
	  cite_please("maragliano06cpl",mtd_data->fplog);
  }
  if(logical.do_pca) cite_please("sutto-2010jctc",mtd_data->fplog);
#ifdef RECONMETAD
  if(reconOn==1) cite_please("tribello-10pnas",mtd_data->fplog);
#endif
  fprintf(mtd_data->fplog,"\n"); 

  disclaimer(mtd_data->fplog);
// flushing output
  fflush(mtd_data->fplog);

}

//-----------------------------------------------------------------------------------------------------------------

void PREFIX read_defaults()
{
  int icv;
 
  colvar.nt_print 		= 10;
  colvar.nconst 		= 0;
  logical.restart_hills 	= 0;
  logical.append                = 0;
  logical.restart_abmd	 	= 0;
  logical.remd 			= 0;
  logical.hrex 			= 0;
  colvar.hrex_energy		= 0.0;
  logical.rpxm			= 0;
  logical.do_hills 		= 0;
  logical.commit 		= 0;
  logical.print 		= 0;
  logical.widthadapt            = 0;
  logical.welltemp              = 0;
  // JFD>
  logical.dicksonian_tempering  = 0;
  logical.transition_tempering  = 0;
  logical.mcgdp_hills      = 0;
  // <JFD
  logical.tamd                  = 0;
  logical.debug                 = 0;
  logical.parallel_hills        = 0;
  logical.norescale             = 0;
#if defined(PLUMED_GROMACS) || defined(DL_POLY) || defined (AMBER)
  logical.parallel_hills        = 1;
#endif
  logical.debug_derivatives     = 0;
  logical.enable_untested_features = 0;
  logical.do_grid               = 0;
  logical.read_grid             = 0;
  logical.write_grid            = 0;
  logical.donot_spline          = 0;
  logical.debug_grid            = 0;
  logical.do_walkers            = 0;
  logical.puckering             = 0;
  logical.path                  = 0;
  logical.energy                = 0;
  logical.read_old_bf           = 0;
  logical.do_external           = 0;
  logical.do_alphabetarmsd      = 0;
  logical.do_sprint             = 0;
  logical.do_steerplan          = 0;
  logical.do_constraint         = 0;
  logical.do_dafed		= 0;   // #### d-AFED
  dafed_control.restart	= 0;   // #### d-AFED checkpointing
  dafed_control.write_freq	= -1;   // #### d-AFED checkpointing
  dafed_control.do_cpt	= 0;   // #### d-AFED checkpointing
  dafed_control.n_respa	= 1;   // #### d-AFED RESPA (default for all non-dafed applications)
  logical.do_inversion          = 0;
  logical.do_pca                = 0;
  sprintf(colvar.hills_label,"\0");
  hills.wwr 			= 0.;
  hills.rate			= 0.;
  hills.max_height              = 0.;
  hills.max_stride              = 0;
  hills.n_hills			= 0;
  hills.nt_hills                = 999999999; 
  hills.nr_hills                = 1; 
  hills.read                    = 0;
  hills.idwalker                = 0;
  hills.nwalkers                = 1;
  nsz                           = 0;
  hills.first_read              = 1;
  colvar.ptmetad_neighbours     = 0;
  colvar.ptmetad_sigma          = 0.0;
  colvar.align_atoms            = 0;
  colvar.align_list             = NULL;
  mtd_data.dump_atoms           = 0;
  mtd_data.dump_list            = 0;
  mtd_data.dump_stride          = 100;
  mtd_data.temp_t		= 0.;
  
  bias_grid.ncv                      = 0;
  bias_grid.nhills                   = 0;
  colvar.pg.list		=NULL;
  colvar.pg.nlist		=0;

  for(icv=0;icv<nconst_max;icv++){
    logical.abmd[icv]	 	= 0;
    logical.steer[icv]          = 0; 
    logical.dafed[icv]		= 0;   // #### d-AFED
    logical.cnstr[icv]          = 0; 
    logical.upper[icv] 		= 0;
    logical.lower[icv] 		= 0;
    logical.interval[icv]       = 0;  
    logical.ureflect[icv]       = 0;
    logical.lreflect[icv]       = 0;
    logical.invert[icv][0]      = 0;
    logical.invert[icv][1]      = 0;
    logical.always[icv]		= 0;
    logical.nlist[icv]          = 0;
    colvar.inv_limit[icv][0]    = -1.e9;
    colvar.inv_limit[icv][1]    = +1.e9;
    colvar.inv_ref[icv]         = 1.6;
    colvar.inv_inv[icv]         = 6;
    colvar.inv_maxww[icv]       = 4;
    cvw.sigma[icv] 		= 0.;
    cvw.upper[icv] 		= 0.;
    cvw.lower[icv] 		= 0.;
    cvw.lsigma[icv] 		= 0.;
    cvw.fwall[icv] 		= 0;
    cvw.uexp[icv] 		= 4;
    cvw.lexp[icv] 		= 4;
    cvw.ueps[icv] 		= 1.;
    cvw.leps[icv] 		= 1.;
    cvw.uoff[icv]               = 0.;
    cvw.loff[icv]               = 0.;
    cvint.lower_limit[icv]        = 0.;   //fahimeh
    cvint.upper_limit[icv]        = 0.;   //fahimeh
    colvar.on[icv] 		= 1;
    colvar.Mss0[icv]            = 0.;
    colvar.M2ss0[icv]           = 0.;
    colvar.type_s[icv]          = 0;
    colvar.logic[icv]           = 0;
    colvar.natoms[icv]          = 0;
    colvar.cell_pbc[icv]        = 0;
    colvar.delta_r[icv]         = -1.; // default synonim of NOHILLS
    colvar.b_treat_independent[icv] = 0;
    bias_grid.min[icv]          = 0.;
    bias_grid.max[icv]          = 0.; 
    bias_grid.lbox[icv]         = 0.;
    bias_grid.dx[icv]           = 0.;
    bias_grid.bin[icv]          = 1;
    bias_grid.minibin[icv]      = 1;
    bias_grid.period[icv]       = 0;
    bias_grid.index[icv]        = 0;
    bias_grid.oldelta[icv]      = 0.;
    // JFD>
    hills.mcgdp_reshape_flag[icv] = 0;
    // <JFD
    //ADW >
    colvar.stoch_sample         = 1.;
    eds.cv_number               = 0;
    // <ADW
    cvsteer.slope[icv]          = 0.;
    cvsteer.annealing[icv]      = 0;
    stopwhen.actmin[icv]        = 0;
    stopwhen.actmax[icv]        = 0;
    steerplan.isactive[icv]        = 0;
  }
  mtd_data.time_offset=0.;
  mtd_data.newcolvarfmt=1;
	// initialize the structure to zero elements
  rmsd_workstruct.maxsize=0;
  rmsd_workstruct.maxsize_secondder=0;
	// initialize the couplingmatrix
	couplingmatrix.is_on=0;

}

//-----------------------------------------------------------------------------------------------------------------

// seek_word WILL BE REMOVED SOON (as soon as it will be replaced everywhere)

int PREFIX seek_word(char **word, const char *wanted)
{
  int i;

  for (i=0;;i++) {
    if (word[i]==NULL) return -1;
    if (strcmp(word[i],wanted)==0) return i;
  }
  return -1;
}

// Added By Paolo to progrssively seek in the input string
int PREFIX seek_word2(char **word, const char *wanted, int is)
{
  int i;

  for (i=is;;i++) {
    if (word[i]==NULL) return -1;
    if (strcmp(word[i],wanted)==0) return i;
  }
  return -1;
}

void PREFIX cite_please (const char* re, FILE *fplog){


 fprintf(fplog, "\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n"); 
 if(!strcmp(re,"laio-parr02pnas")){
    fprintf(fplog, "  A. Laio and M. Parrinello\n");
    fprintf(fplog, "  Escaping free energy minima\n");
    fprintf(fplog, "  Proc. Natl. Acad. Sci. USA. 2002 vol. 99 (20) pp. 12562-6\n");
 } else if(!strcmp(re,"bono+09cpc")){
    fprintf(fplog, "  M. Bonomi, D. Branduardi, G. Bussi, C. Camilloni, D. Provasi, P. Raiteri, \n");
    fprintf(fplog, "  D. Donadio, F. Marinelli, F. Pietrucci, R. A. Broglia and M. Parrinello \n");
    fprintf(fplog, "  PLUMED: a portable plugin for free-energy calculations with molecular dynamics\n");
    fprintf(fplog, "  Comp. Phys. Comm. 2009 vol. 180 (10) pp. 1961-1972 \n");
 } else if(!strcmp(re,"pian-laio07jpcb")){
    fprintf(fplog, "  S. Piana and A. Laio\n");
    fprintf(fplog, "  A Bias-Exchange Approach to Protein Folding \n");
    fprintf(fplog, "  J. Phys. Chem. B. 2007 vol. 111 (17) pp. 4553-9\n");
 } else if(!strcmp(re,"bard+08prl")){
    fprintf(fplog, "  A. Barducci, G. Bussi and M. Parrinello\n");
    fprintf(fplog, "  Well-Tempered Metadynamics: A Smoothly Converging and Tunable Free-Energy Method \n");
    fprintf(fplog, "  Phys. Rev. Lett. 2008 vol. 100 (2) pp. 020603 \n");
 } else if(!strcmp(re,"buss+06jacs")){
    fprintf(fplog, "  G. Bussi, F.L. Gervasio, A. Laio and M. Parrinello \n");
    fprintf(fplog, "  Free-energy landscape for beta hairpin folding from combined parallel tempering and metadynamics\n");
    fprintf(fplog, "  J. Am. Chem. Soc. 2006 vol. 128 (41) pp. 13435-41 \n");
 } else if(!strcmp(re,"bran+07jcp")){
    fprintf(fplog, "  D. Branduardi, F.L. Gervasio and M. Parrinello \n");
    fprintf(fplog, "  From A to B in free energy space\n");
    fprintf(fplog, "  Jour. Chem. Phys. 2007 vol. 126 (5) pp. 054103\n");
 } else if(!strcmp(re,"rait+06jpcb")){
    fprintf(fplog, "  P. Raiteri, A. Laio, F.L. Gervasio, C. Micheletti and M. Parrinello \n");
    fprintf(fplog, "  Efficient Reconstruction of Complex Free Energy Landscapes by Multiple Walkers Metadynamics \n"); 
    fprintf(fplog, "  J. Phys. Chem. B. 2006 vol. 110 (8) pp. 3533-3539 \n");
 } else if(!strcmp(re,"sega+09jcp")){
    fprintf(fplog, "  M. Sega, E. Autieri and F. Pederiva\n");
    fprintf(fplog, "  On the Calculation of Puckering Free Energy Surfaces \n"); 
    fprintf(fplog, "  J. Chem. Phys. 2009 vol. 130 (22) pp. 225102 \n");
 } else if(!strcmp(re,"pietrucci+09jctc")){
    fprintf(fplog, "  F. Pietrucci and A. Laio\n");
    fprintf(fplog, "  A collective variable for the efficient exploration of protein beta-structures with metadynamics: application to SH3 and GB1\n");
    fprintf(fplog, "  J. Chem. Theory Comput. 2009 vol. 5(9) pp. 2197 \n");
 } else if(!strcmp(re,"pietrucci+11prl")){
    fprintf(fplog, "  F. Pietrucci and W. Andreoni\n");
    fprintf(fplog, "  Graph theory meets ab initio molecular dynamics: atomic structures and transformations at the nanoscale\n");
    fprintf(fplog, "  Phys. Rev. Lett. 2011 vol. 107(8) pp. 085504\n");
 } else if(!strcmp(re,"marinell-crespo10")){
   fprintf(fplog, "  Y. Crespo, F. Marinelli, F. Pietrucci, A. Laio\n");
   fprintf(fplog, "  Metadynamics convergence law in a multidimensional system\n");
   fprintf(fplog, "  Phys. Rev. E 2010 vol. 81(5) pp. 055701 \n");
 } else if(!strcmp(re,"tribello-10pnas")){
   fprintf(fplog, "  G. A. Tribello, M. Ceriotti and M. Parrinello\n");
   fprintf(fplog, "  A self-learning algorithm for based molecular dynamics\n");
   fprintf(fplog, "  Proc. Natl. Acad. Sci. U.S.A. 2010 vol 107(41) pp. 17509-17514\n");  
 } else if(!strcmp(re,"sutto-2010jctc")){
   fprintf(fplog, "  L. Sutto, M. D'Abramo and F.L. Gervasio\n");
   fprintf(fplog, "  Comparing the Efficiency of Biased and Unbiased Molecular Dynamics\n");
   fprintf(fplog, "  in Reconstructing the Free Energy Landscape of Met-Enkephalin\n");
   fprintf(fplog, "  J. Chem. Theory Comput. 2010 vol 6(12) pp.3640-3646\n");  
 } else if(!strcmp(re,"abrams08jpcb")){         // #### d-AFED
    fprintf(fplog, "  J. B. Abrams and M. E. Tuckerman\n");
    fprintf(fplog, "  Efficient and Direct Generation of Multidimentional Free Energy Surfaces via Adiabatic Dynamics without Coordinate Transformations \n");
    fprintf(fplog, "  J. Phys. Chem. B 2008 vol. 112 pp. 15742-15757 \n");
 } else if(!strcmp(re,"maragliano06cpl")){      // #### d-AFED
    fprintf(fplog, "  L. Maragliano and E. Vanden-Eijnden\n");
    fprintf(fplog, "  A Temperature Accelerated Method for Sampling Free Energy and Determining Reaction Pathways in Rare Events Simulations \n");
    fprintf(fplog, "  Chem. Phys. Lett. 2006 vol. 426 pp. 168-175 \n");
 // JFD>
 } else if(!strcmp(re,"dickson2011pre")){
    fprintf(fplog, "  B. M. Dickson\n");
    fprintf(fplog, "  Approaching a parameter-free metadynamics \n");
    fprintf(fplog, "  Phys. Rev. E 2011 vol. 84 pp. 037701 \n");
 } else if(!strcmp(re,"mcgovern2013jcp")){
    fprintf(fplog, "  M. McGovern and J. J. De Pablo\n");
    fprintf(fplog, "  A boundary correction algorithm for metadynamics in multiple dimensions.\n");
    fprintf(fplog, " J. Chem. Phys. 2013 vol. 139 pp. 084102\n");
 // <JFD
 } else {
    assert(1); // wrong bib name
 }

 fprintf(fplog, "-------- -------- --- Thank You --- -------- --------\n\n"); 
};

void PREFIX disclaimer (FILE *fplog){

 fprintf(fplog,"** PLUMED is free software: you can redistribute it and/or modify \n");
 fprintf(fplog,"** it under the terms of the GNU Lesser General Public License as published by \n");
 fprintf(fplog,"** the Free Software Foundation, either version 3 of the License, or \n");
 fprintf(fplog,"** (at your option) any later version. \n\n");
 fprintf(fplog,"** PLUMED is distributed in the hope that it will be useful,\n");
 fprintf(fplog,"** but WITHOUT ANY WARRANTY; without even the implied warranty of \n");
 fprintf(fplog,"** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n");
 fprintf(fplog,"** GNU Lesser General Public License for more details. \n\n");
 fprintf(fplog,"** You should have received a copy of the GNU Lesser General Public License\n");
 fprintf(fplog,"** along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.  \n\n");
 fprintf(fplog,"** For more info, see:  http://www.plumed-code.org \n");
 fprintf(fplog,"** or subscribe to plumed-users@googlegroups.com \n\n");
#ifdef NAMD
 fprintf(fplog,"                              WARNING!! \n");
 fprintf(fplog,"** Starting from version 2.7, NAMD has its own module for collective  \n");
 fprintf(fplog,"** variable-based calculations including metadynamics, adaptive biasing \n");
 fprintf(fplog,"** force method, umbrella sampling and steered molecular dynamics. \n");
 fprintf(fplog,"** Please, have a look at the NAMD manual for more info. \n\n");
#endif
}; 


//.................................
//.. HERE WE HAVE THE NEW PARSER ..
//.................................

// very long lines allowed
#define PLUMED_LINEMAX  50000

// split a line into words
int PREFIX plumed_get_words(char* line,char*** words){
  char* ww;
  int i;
  (*words)=NULL;
  
  ww=strtok(line," \t\n"); if(ww==NULL) return 0;
  srenew(*words,1);
  (*words)[0]=ww;
  for(i=1;(ww=strtok(NULL," \t\n"));i++){
    srenew(*words,i+1);
    (*words)[i]=ww;
  };
  return i;
};

void PREFIX plumed_error(const char*s){
  fprintf(stderr,"!!!!! PLUMED ERROR: %s\n",s);
  fprintf(stderr,"!!!!! ABORTING RUN \n");
  if(mtd_data.fplog) {
    fprintf(mtd_data.fplog,"!!!!! PLUMED ERROR: %s\n",s);
    fprintf(mtd_data.fplog,"!!!!! ABORTING RUN \n");
    fflush(mtd_data.fplog);
  };
  EXIT();
};

void PREFIX plumed_warn(const char*s){
  if(mtd_data.fplog) {
    fprintf(mtd_data.fplog,"|- : %s\n",s);
    fflush(mtd_data.fplog);
  };
};

// parse a word:
//   if the word ends with postfix, delete the postfix and return 1
//   otherwise return 0
// example:
//   char* word; word=malloc(100); strcpy(word,"pippo->");
//   plumed_parse_word(word,"<-"); // returns 0 and does not change word
//   plumed_parse_word(word,"->"); // returns 1 and changes word to "pippo"
int plumed_parse_word(char* word,const char* postfix){
  int lword;
  int lpostfix;
  lword=strlen(word);
  lpostfix=strlen(postfix);
  if(lword<lpostfix) return 0;
  if(strcmp(& word[lword-lpostfix],postfix)) return 0;
  word[lword-lpostfix]=0;
  return 1;
};

int PREFIX plumed_atoi(const char* word){
  int n,i;
  char* buf;
  n=-1;
  sscanf(word,"%i%n",&i,&n);
  if(n!=strlen(word)) {
    snew(buf,strlen(word)+100);
    sprintf(buf,"parsing integer %s\n",word);
    plumed_error(buf);
  }
  return i;
};

double PREFIX plumed_atof(const char* word)
{
 double x;
 int n, found=0; 

 n=-1;
 sscanf(word,"%lf%n",&x,&n);
 if(n!=strlen(word)){
 if(strcmp(word,"pi")==0)    { x=M_PI;      found=1; } 
 if(strcmp(word,"+pi")==0)   { x=M_PI;      found=1; } 
 if(strcmp(word,"-pi")==0)   { x=-1.*M_PI;  found=1; }
 if(strcmp(word,"2pi")==0)   { x=M_2PI;     found=1; }
 if(strcmp(word,"+2pi")==0)  { x=M_2PI;     found=1; }
 if(strcmp(word,"-2pi")==0)  { x=-1.*M_2PI; found=1; }
 if(found==0) plumed_error("Special symbol not recognized. Please, read the manual for accepted symbols.");
 }
 return x;
};
/*
 get group takes a word  and looks if it's in <mygroup> format:
 if it's so it reads the parsed input and looks for a group specified as

 mygroup->       
   4567 678 678 567 4567
   67 89 678 34 234
   5 7
 mygroup-<       

  IT returns the number of found field for the group ( in the case above 12 ) 
  it places the member of the group in  a vector vec , starting from position n
  ( so the new positions will be stored in vec[i] i=n,i<n+j,i++  ) and reallocate the 
  vector if necessary    

  if the word is not in <mygroup> format, it interprets it as a single number and adds it to the atoms list
   
*/

int PREFIX plumed_get_group(const char *word,int **atoms,int n,t_plumed_input* input,FILE *log){
  int lword,justoneatom,foundgroup,nadd;
  int* toadd;
  lword=strlen(word);
  foundgroup=0;
// check for group syntax
  if(lword>2) if(word[0]=='<' && word[lword-1]=='>') foundgroup=1;
// if not group, just read the atomx index;
  if(!foundgroup){
    justoneatom=plumed_atoi(word)-1;
    nadd=1;
    toadd=&justoneatom;
  } else {
// if group, search for it on the list
    int igroup;
    int found;
    char* groupname;
    snew(groupname,strlen(word)-1);
    strncpy(groupname,& word[1],strlen(word)-2);
    found=0;
    for(igroup=0;igroup<input->ngroups;igroup++){
      if(!strcmp(groupname,input->groupnames[igroup])){
        found=1;
        break;
      }
    }
    sfree(groupname);
    if(!found) plumed_error("group not found");
    nadd=input->natoms[igroup];
    toadd=input->atoms[igroup];
  }
  srenew((*atoms),n+nadd);
  int i;
  for(i=0;i<nadd;i++) (*atoms)[i+n]=toadd[i];
  return nadd;
  
};

// routine to parse input file
// * remove comments
// * join lines with continuation
// * find and stores the group definitions
// * save everything which is not a group in the array input->words[iline][iword]
//   where iline=0...(input->nlines-1) and iword=0...(input=->nwords[iline)
// * line numbers are preserved to allow a better error reporting
void PREFIX plumed_read_input(t_plumed_input* input,FILE* file,FILE* log){
  char* line;
  int i;
  int iline,iword;
  char** words;
  int nwords;
  char* inside_group;
  int inside_loop,loop_start,loop_end,loop_stride;

// initial values
  input->nlines=0;
  input->nwords=NULL;
  input->words=NULL;
  input->ngroups=0;
  input->groupnames=NULL;
  input->natoms=NULL;
  input->atoms=NULL;

  inside_group=NULL;
  inside_loop=0;
  loop_start=0;
  loop_end=0;

  snew(line,PLUMED_LINEMAX);

  while(fgets(line,PLUMED_LINEMAX,file)){

   iline=input->nlines;
   input->nlines++;
   srenew(input->nwords,input->nlines);
   input->nwords[iline]=0;
   srenew(input->words,input->nlines);
   input->words[iline]=NULL;

// merge lines ending with "backslash" or "ampersand"
   int linelength;
   linelength=strlen(line);
   if(linelength>1) while(line[linelength-2]=='\\' || line[linelength-2]=='&'){
     if(!fgets(&line[linelength-2],PLUMED_LINEMAX-linelength+2,file))
       plumed_error("last line is not ending");
     linelength=strlen(line);
// append an empty line
// in this way the line count corresponds to the file (better for error reporting)
     iline=input->nlines;
     input->nlines++;
     srenew(input->nwords,input->nlines);
     input->nwords[iline]=0;
     srenew(input->words,input->nlines);
     input->words[iline]=NULL;
// AN EXTRA EMPTY WORD IS ADDED TO BE COMPATIBLE WITH OLDER seek_word
       srenew(input->words[iline],1);
       input->words[iline][0]=NULL;
   }

// Remove comments (beginning with sharp or esclamation)
    for(i=0;line[i];i++) if(line[i]=='#' || line[i]=='!') line[i]=0;

// Split into words:
    nwords=plumed_get_words(line,&words);

// Check for ENDMETA or ENDPLUMED
    if(nwords>0) if(!strcmp(words[0],"ENDMETA") || !strcmp(words[0],"ENDPLUMED")) {
      free(words);
      break;
    }

// loop over all the input words
    for(iword=0;iword<nwords;iword++){

// begin group
      if(plumed_parse_word(words[iword],"->")){
        int igroup;
        if(inside_group) plumed_error("nested groups are not allowed");
        if(iword>0) plumed_error("a group cannot begin in the middle of a line");
        igroup=input->ngroups;
        input->ngroups++;
//   store group name
        srenew(input->groupnames,input->ngroups);
        snew(input->groupnames[igroup],strlen(words[iword])+1);
        strcpy(input->groupnames[igroup],words[iword]);
//   initialize its atom list
        srenew(input->natoms,input->ngroups);
        input->natoms[igroup]=0;
        srenew(input->atoms,input->ngroups);
        input->atoms[igroup]=NULL;
        inside_group=input->groupnames[igroup];

        fprintf(log,"|- GROUP FOUND: %s\n",inside_group);

//   check if other groups with the same name have been defined
        for(i=0;i<igroup;i++) if(!strcmp(inside_group,input->groupnames[i]))
          plumed_error("two groups cannot have the same name");

// end group
      } else if(plumed_parse_word(words[iword],"<-")){
        int igroup,iatom;
        igroup=input->ngroups-1;
        if(!inside_group) plumed_error("end group without begin group");
        if(strcmp(inside_group,words[iword])) plumed_error("end group different from begin group");
        if(inside_loop==1 || inside_loop==2) plumed_error("wrong LOOP syntax");
// this is for backward compatibility with "LOOP 1 10", without stride
        if(inside_loop==3){
          for(i=loop_start;i<=loop_end;i++){
            iatom=input->natoms[igroup];
            input->natoms[igroup]++;
            srenew(input->atoms[igroup],input->natoms[igroup]);
            input->atoms[igroup][iatom]=i-1;
          }
        };

// log the list of members
        fprintf(log,"|- GROUP MEMBERS: ");
        for(i=0;i<input->natoms[igroup];i++){
          if((i+1)%20==0) fprintf(log,"\n|-                ");
          fprintf(log," %i",input->atoms[igroup][i]+1);
        }
        fprintf(log,"\n");
        inside_group=NULL;

// if we are within a group definition, add atom of check for loop syntax
      } else if(inside_group) {
        int i,igroup,iatom;
        igroup=input->ngroups-1;
// NOTE: this should be triggered only by a LOOP keyword inside a group definition
//       it should allow for a hypothetical LOOP keyword in a standard directive
        if(!strcmp("LOOP",words[iword])){
          inside_loop=1;
        } else if(inside_loop==1){
          loop_start=plumed_atoi(words[iword]);
          inside_loop=2;
        } else if(inside_loop==2){
          loop_end=plumed_atoi(words[iword]);
          inside_loop=3;
        } else if(inside_loop==3){
          loop_stride=plumed_atoi(words[iword]);
          for(i=loop_start;i<=loop_end;i+=loop_stride){
            iatom=input->natoms[igroup];
            input->natoms[igroup]++;
            srenew(input->atoms[igroup],input->natoms[igroup]);
            input->atoms[igroup][iatom]=i-1;
          }
          inside_loop=0;
        } else {
//   add a single atom to the list
          i=plumed_atoi(words[iword]);
          iatom=input->natoms[igroup];
          input->natoms[igroup]++;
          srenew(input->atoms[igroup],input->natoms[igroup]);
          input->atoms[igroup][iatom]=i-1;
        }
// if we are on a normal line, just copy the word
      } else {
        int iw;
        iw=input->nwords[iline];
        input->nwords[iline]++;
        srenew(input->words[iline],input->nwords[iline]);
// // AN EXTRA EMPTY WORD IS ADDED TO BE COMPATIBLE WITH OLDER seek_word
          srenew(input->words[iline],input->nwords[iline]+1);
          input->words[iline][iw+1]=NULL;
        snew(input->words[iline][iw],strlen(words[iword])+1);
        strcpy(input->words[iline][iw],words[iword]);
      };
    }

// Finally delete word pointer for this line
    free(words);
  };

// This buffer is not needed anymore
  sfree(line);
  
// DEBUG
//  for(iline=0;iline<input->nlines;iline++) for(iword=0;iword<input->nwords[iline];iword++)
//  fprintf(log,"%i %i : '%s'\n",iline,iword,input->words[iline][iword]);
};

// Deallocate memory
void PREFIX plumed_clear_input(t_plumed_input*input){
  int i,j;
  for(i=0;i<input->nlines;i++) for(j=0;j<input->nwords[i];j++) sfree(input->words[i][j]);
  sfree(input->nwords);
  for(i=0;i<input->nlines;i++) sfree(input->words[i]);
  sfree(input->words);
  for(i=0;i<input->ngroups;i++) sfree(input->groupnames[i]);
  for(i=0;i<input->ngroups;i++) sfree(input->atoms[i]);
  sfree(input->groupnames);
  sfree(input->atoms);
  sfree(input->natoms);
  input->nlines=0;
  input->nwords=NULL;
  input->words=NULL;
  input->ngroups=0;
  input->groupnames=NULL;
  input->natoms=NULL;
  input->atoms=NULL;
}

void PREFIX couple2list( struct coupling_ll **first_elem ,int *at1,int nat1,int *at2,int nat2){
     int i,j,k;
     struct coupling_ll *newelem,mycouple; 
     newelem= (struct coupling_ll *)malloc(sizeof( struct coupling_ll)); 

     printf("SUMMARY************************** %p \n",first_elem);
     // this add the element to the linked list
     (* newelem).nat1=nat1;   
     (* newelem).at1=(int *)malloc(nat1*sizeof(int));   
     for(i=0;i<nat1;i++){ (* newelem).at1[i]=at1[i];}
     (* newelem).nat2=nat2;   
     (* newelem).at2=(int *)malloc(nat2*sizeof(int));   
     for(i=0;i<nat2;i++){ (* newelem).at2[i]=at2[i];}

     // the address of newelem is pointing to first elem
     (* newelem).next_elem= (*first_elem);
     // now the pointer first elem is pointing to the new elem  
     (* first_elem) = newelem; 
     printf("AT1 "); 
     for(i=0;i<(**first_elem).nat1;i++){
          printf(" %d ",(**first_elem).at1[i]);
     } 
     printf("\n");
     printf("AT2 "); 
     for(i=0;i<(**first_elem).nat2;i++){
         printf(" %d ",(**first_elem).at2[i]);
     } 
     printf("\n");
     printf("ENDSUMMARY************************** %p \n",first_elem);
};
void PREFIX freecouple2list( struct coupling_ll **first_elem ){
     struct coupling_ll *newelem; 
        // this add the element to the linked list
         while( (* first_elem)!=NULL){
             free((**first_elem).at1);
             free((**first_elem).at2);
             newelem= (* first_elem);
             (* first_elem)=(* newelem).next_elem;
             free(newelem);
         }
};
void PREFIX scancouple( struct coupling_ll *first_elem ){
     int i,j,k;
     struct coupling_ll *ptr; 
     ptr=first_elem;
     printf("SCANCOUPLE*************************\n");
     while(ptr!=NULL){
        printf("NEWCOUPLE************************* %p\n",ptr);
        printf("NAT1 %d ",(*ptr).nat1);  
        //EXIT(); 
        for(j=0;j<(*ptr).nat1;j++){
           k=(*ptr).at1[j];
           printf(" AT %d ",k);
        }
        printf("\n") ;
        if((*ptr).nat2){
            printf("NAT2 %d ",(*ptr).nat2);  
            for(j=0;j<(*ptr).nat2;j++){
               k=(*ptr).at2[j];
                printf(" AT %d ",k);
            }
        }
        printf("\n");
        ptr=ptr->next_elem;
     }   
     printf("ENDSCANCOUPLE*************************\n");
 
};
void PREFIX setup_projections(struct proj_grad_s *proj ){
    int i,j,k,dimension,ncv;
    int ii,jj,kk,ll,iii,jjj;
    int nat1,*at1;
    int nat2,*at2;
    int *skip; 
    ncv=proj->nlist;
    dimension=mtd_data.natoms; // total number of atoms
    skip=(int *)malloc(dimension*sizeof(int)); 
    at1=(int *)malloc(dimension*sizeof(int)); 
    at2=(int *)malloc(dimension*sizeof(int)); 

    proj->matrix=(real **)malloc(colvar.nconst*sizeof(real *));
    for(i=0;i<colvar.nconst;i++){
       proj->matrix[i]=(real *)malloc(colvar.nconst*sizeof(real));
    }
 
    fprintf(mtd_data.fplog,"|-PROJ_GRAD: TOTAL DIMENSION %d\n",dimension);
    proj->couple=(struct el_couple *)malloc((ncv*(ncv-1)/2)*sizeof( struct el_couple )); //one element for couple         
    for(i=0;i<(ncv*(ncv-1)/2);i++)(proj->couple[i]).first_elem=NULL; // set each pointer for the linked list  to null 
    k=0;// progressive for couple counting
    for(iii=0;iii<proj->nlist-1;iii++){
      i=proj->list[iii]; 
      for(jjj=iii+1;jjj<proj->nlist;jjj++){
          j=proj->list[jjj];
          proj->couple[k].cv1=i;
          proj->couple[k].cv2=j;
          for(ii=0;ii<dimension;ii++){skip[ii]=0;} 
          for(ii=0;ii<colvar.natoms[i];ii++){
               nat1=0;nat2=0;
               kk=colvar.cvatoms[i][ii]; // get the index of atom 
               if(skip[ii]){
                   fprintf(mtd_data.fplog,"|-PROJ_GRAD: SKIPPING ATOM %d \n",ii); 
               }else{  
                   // find within the  other set  
                   for(jj=0;jj<colvar.natoms[j];jj++){ // find within the other cv
                           ll=colvar.cvatoms[j][jj]; 
          //                 printf("ATOM1 %d ATOM2 %d \n",kk,ll);
                           if(kk==ll){// found common indexes
                              at2[nat2]=jj;            
                              nat2++;
          //                    printf("ATOM1=ATOM2 \n");
                           }
                   } 
                   if(nat2){ // makes sense only when there is at least one nat2
                        // find within the  same set ( for additive cv ) and exclude computation 
                        for(jj=ii;jj<colvar.natoms[i];jj++){ // find within the other cv
                            ll=colvar.cvatoms[i][jj]; 
                            if(kk==ll){// found common indexes
                               at1[nat1]=jj;            
                               nat1++;
                               skip[jj]=1;
                            }
                        }
                        // transfer it to the right vectors of the linked list 
                   }  
                   if(nat2){
			couple2list( &((proj->couple[k]).first_elem) ,at1,nat1,at2,nat2);
		   }	
               }
          } 
          // freecouple2list(&((proj.couple[k]).first_elem) );
          scancouple((proj->couple[k]).first_elem);
          // increment couple counter

          k++;
      }  
    }

     proj->ncouples=k;
    // printf("NCOUPLES FOUND %d\n",proj.ncouples);
     // DIAGONAL CONTRIBUTION
     fprintf(mtd_data.fplog,"|-PROJ_GRAD: DIAGONAL CONTRIBUTION \n"); 
     proj->diagonal=(struct el_diagonal *)malloc(ncv*sizeof( struct el_diagonal )); //one element for couple         
     for(i=0;i<ncv;i++)(proj->diagonal[i]).first_elem=NULL; // set each pointer for the linked list  to null 
     for(i=0;i<ncv;i++){
               for(ii=0;ii<dimension;ii++)skip[ii]=0; 
               // look for all the atoms which  have the same index
               for(ii=0;ii<colvar.natoms[i];ii++){
                   nat1=0;nat2=0;
                   kk=colvar.cvatoms[i][ii]; // get the index of atom 
                   if(skip[ii]){
                       fprintf(mtd_data.fplog,"|-PROJ_GRAD: SKIPPING ATOM %d\n",ii); 
                   }else{  
                   // find within the  same set ( for additive cv ) and exclude computation 
                      for(jj=ii;jj<colvar.natoms[i];jj++){ // find within the other cv
                         ll=colvar.cvatoms[i][jj]; 
                         if(kk==ll){// found common indexes
                            at1[nat1]=jj;            
                            nat1++;
                            skip[jj]=1;
                         }
                      }
                   }
                   // transfer it to the right vectors of the linked list 
                   if(nat1)couple2list( &((proj->diagonal[i]).first_elem) ,at1,nat1,at2,nat2);
               }
        //       freecouple2list(&((proj.diagonal[i]).first_elem) );
               scancouple((proj->diagonal[i]).first_elem);
     }
    free(skip);
    free(at1);
    free(at2);
    // now create the header for COLVAR file
    sprintf(proj->log," XX "); 
    char cvlog[100];
    int cv1,cv2; 
    ncv=proj->nlist;
    for(i=0;i<ncv;i++){
        cv1=proj->list[i];
        for(j=i;j<ncv;j++){
           cv2=proj->list[j];
           sprintf(cvlog," NABLACV%d_DOT_NABLACV%d ",cv1+1,cv2+1);
 	   strcat(proj->log,cvlog); 
        } 
    }
    return; 
}
void PREFIX calc_projections(struct proj_grad_s  *proj ){
    int ncv,cv1,cv2,i,j,k,ii,jj;
    real scal,grad;
    real tmp1x,tmp1y,tmp1z;
    real tmp2x,tmp2y,tmp2z;
    struct coupling_ll *ptr;
    ncv= proj->nlist;
    for(i=0;i<colvar.nconst;i++){
        for(j=i;j<colvar.nconst;j++){
           (proj->matrix[i][j])=0.;
        }
    }
    for(i=0;i<proj->ncouples;i++){// loop over all the couples
            //printf("COUPLE %d ADDR %d\n",i,proj->couple[i].first_elem);
            cv1=proj->couple[i].cv1;
            cv2=proj->couple[i].cv2;
            scal=0.;
            // use the linked list to calculate grad(cv1) dot grad(cv2)
            ptr=proj->couple[i].first_elem;
            while(ptr!=NULL){
               tmp1x=0.;tmp1y=0.;tmp1z=0.; 
               for(j=0;j<(*ptr).nat1;j++){
                  k=(*ptr).at1[j];
                  //printf("KK %d \n",k);
                  tmp1x+=colvar.myder[cv1][k][0];
                  tmp1y+=colvar.myder[cv1][k][1];
                  tmp1z+=colvar.myder[cv1][k][2];
               }
               tmp2x=0.;tmp2y=0.;tmp2z=0.; 
               for(j=0;j<(*ptr).nat2;j++){
                  k=(*ptr).at2[j];
                  //printf("MM %d \n",k);
                  tmp2x+=colvar.myder[cv2][k][0];
                  tmp2y+=colvar.myder[cv2][k][1];
                  tmp2z+=colvar.myder[cv2][k][2];
               }
               scal+=tmp1x*tmp2x+
                     tmp1y*tmp2y+
                     tmp1z*tmp2z;
               ptr=ptr->next_elem;
            }   
            proj->matrix[cv1][cv2]=scal; 
         //   printf("II %d JJ %d VV %f \n",cv1,cv2,scal);
    }
    for(i=0;i<ncv;i++){
               cv1=proj->list[i];
               grad=0.0;
               ptr=proj->diagonal[i].first_elem;
               while(ptr!=NULL){
                   tmp1x=0.;tmp1y=0.;tmp1z=0.; 
                   for(j=0;j<(*ptr).nat1;j++){
                      k=(*ptr).at1[j];
                      //printf("NN %d \n",k);
                      tmp1x+=colvar.myder[cv1][k][0];
                      tmp1y+=colvar.myder[cv1][k][1];
                      tmp1z+=colvar.myder[cv1][k][2];
                   }
                   grad+=tmp1x*tmp1x;
                   grad+=tmp1y*tmp1y;
                   grad+=tmp1z*tmp1z;
                   ptr=ptr->next_elem;
               }
               proj->matrix[cv1][cv1]=grad; 
        //    printf("II %d JJ %d VV %f \n",cv1,cv1,grad);
     }
     for(i=0;i<ncv-1;i++){
        ii=proj->list[i];
        for(j=i+1;j<ncv;j++){
           jj=proj->list[j];
           proj->matrix[jj][ii]=proj->matrix[ii][jj];
        }
     }
 //    printf("MMM \n");
//        printf("MMM ");
     sprintf(proj->log," PROJ_GRAD "); 
     char cvlog[100];
     for(i=0;i<ncv;i++){
        ii=proj->list[i];
        for(j=i;j<ncv;j++){
           jj=proj->list[j];
           sprintf(cvlog," %12.6f ",proj->matrix[ii][jj]); 
           //printf(" %f ",proj->matrix[ii][jj]);
           strcat(proj->log,cvlog);
        }
     }
     //  printf("\n");
    return;
}

int PREFIX read_steerplan(char **word,int count,t_plumed_input *input,int *iline,FILE *fplog) {
    int i, iw, iw2, iat,j,cvs[100];
    char string[400],chr[3];
    int help;
    help=0;
    fprintf(fplog,"|-READING STEERPLAN\n");  
    iw2 = seek_word(word,"STEERPLAN");

    if (iw2 >=0){
       load_steerplan(word[iw2+1],&steerplan,fplog);             
    } else{ fprintf(fplog,"|- NEEDED FILE KEYWORD FOR STEERPLAN\n"); help=1;}

    fprintf(fplog,"|-END READING STEERPLAN\n");  
    return 1;
};

int PREFIX  load_steerplan  (char *filename,struct steerplan_s *steerplan, FILE *fplog){
    int i,j,k,l,m,jj,ii,kk,nw,ncvs,*maxcvs,endcv;
    FILE *fp;  
    char string[200],*str,**words,string2[200],string3[200];
    snew(maxcvs,100);
    for(i=0;i<100;i++)maxcvs[i]=-1; // default: doesn't have any cv
    fprintf(fplog,"|--LOADING STEERPLAN\n");  
    // open and verify the file
    fp=fopen(filename,"r");
    if (fp == NULL) {
       char buf[1024];
       sprintf(buf,"UNSUCCESSFULL OPENING FOR FILE %s",filename);
       plumed_error(buf);  
    }
    j=0;k=0;
    while(1){
        newline1:
        str=fgets(string,200,fp);
        //first word marks a comment? skip it... 
        nw=plumed_get_words(string,&words); 
        if(!strcmp("#",words[0])) goto newline1; 
        if(nw==0) goto newline1;
        // find  all the occurrency of CV and define the maximum of cvs involved in the 
        // steerplan
        for(i=0;i<nw;i++){
           if(!strcmp("CV",words[i]) && ( maxcvs[atoi(words[i+1])-1]<0 )  ){ maxcvs[atoi(words[i+1])-1]=k;k++;  }; 
        } 

        if(str==NULL)break;
        if (feof(fp))break;
        j++;
    }
    // now k has the number of actions , j has the number of cvs: allocation
    fclose(fp);
    // allocate the steerplan
    steerplan->ncvs=k;
    // create an  index for used cvs
    fprintf(mtd_data.fplog,"|--found %d MAXIMUM CVS INVOLVED IN THE STEERPLAN\n",steerplan->ncvs); 
    (steerplan->totstages)=j;
    fprintf(mtd_data.fplog,"|--found %d STEERACTIONS\n",steerplan->totstages); 
    snew(steerplan->actions,steerplan->totstages); 
    // allocate the for each action allocate the number of cvs 
    for(i=0;i<steerplan->totstages;i++){
       snew(steerplan->actions[i].activecv,steerplan->ncvs);  
    }
    snew(steerplan->actualcv,steerplan->ncvs);
    // clean to defaults
    for(i=0;i<steerplan->totstages;i++){
        for(j=0;j< steerplan->ncvs;j++){
                   steerplan->actions[i].activecv[j].k=-1; // signal for unassigned action in this section  
                   steerplan->actions[i].activecv[j].type=1; 
                   for(k=0;k<100;k++){ if(maxcvs[k]==j){ steerplan->isactive[k]=1;steerplan->actions[i].activecv[j].ncv=k; break;}} 
                   steerplan->actions[i].activecv[j].wildcardpos=0; 
                   steerplan->actions[i].activecv[j].wildcardk=0; 
        }
    }
    // load the steerplan 
    jj=0;
    fp=fopen(filename,"r");
    while(1){
        newline:
        str=fgets(string,200,fp);
        if(str==NULL)break;
        if (feof(fp))break;
        i=plumed_get_words(string,&words); 
        if(!strcmp("#",words[0])) goto newline; 
        if(i==0)goto newline; 
//
// format :
//  time   CV 1  k pos type CV 2 k pos type CV 3 k pos type
//
            steerplan->actions[jj].t=atof(words[0]);
            // clean to default 
           // split the number of needed cvs
            for(j=1;j<i-1;j++){
              if(!strcmp(words[j],"CV")){
                   endcv=i; // default is the end of the line
                   for(k=j+1;k<i;k++){
                        if(!strcmp(words[k],"CV")){endcv=k; break;};
                   }
                   //fprintf(mtd_data.fplog,"|--found atomic from %d to %d \n",j,endcv); 
                   //parse only this cv directive (too atomized?? )
                   if(endcv-j==5){
                       l=atoi(words[j+1])-1;// index of the cv 
                       // find it in the array  
                       m=maxcvs[l];//  fprintf(mtd_data.fplog,"|- collocating action %d  cv %d into array pos %d \n",jj,l+1,m); 
                       steerplan->actions[jj].activecv[m].ncv=l;
                       if(!strcmp(words[j+2],"*")){
                          steerplan->actions[jj].activecv[m].wildcardk=1;
                       } else {
                          steerplan->actions[jj].activecv[m].k=atof(words[j+2]);
                       } 
                       if(!strcmp(words[j+3],"*")){
                          steerplan->actions[jj].activecv[m].wildcardpos=1;
                       } else {
                          steerplan->actions[jj].activecv[m].pos=atof(words[j+3]);
                       } 
                       if(!strcmp(words[j+4],"CENTRAL")){
                          steerplan->actions[jj].activecv[m].type=1;
                       }else if(!strcmp(words[j+4],"POSITIVE")){
                          steerplan->actions[jj].activecv[m].type=2;
                       }else if(!strcmp(words[j+4],"NEGATIVE")){
                          steerplan->actions[jj].activecv[m].type=3;
                       }else {
                           plumed_error("The only types of constraints accepted in the steerplan are CENTRAL POSITIVE NEGATIVE (default is CENTRAL of course)");  
                       }
                   }else if(endcv-j==4){
                       l=atoi(words[j+1])-1;// index of the cv 
                       // find it in the array  
                       m=maxcvs[l]; // fprintf(mtd_data.fplog,"|- collocating action %d  cv %d into array pos %d \n",jj,l+1,m); 
                       steerplan->actions[jj].activecv[m].ncv=l;
                       if(!strcmp(words[j+2],"*")){
                          steerplan->actions[jj].activecv[m].wildcardk=1;
                       } else {
                          steerplan->actions[jj].activecv[m].k=atof(words[j+2]);
                       } 
                       if(!strcmp(words[j+3],"*")){
                          steerplan->actions[jj].activecv[m].wildcardpos=1;
                       } else {
                          steerplan->actions[jj].activecv[m].pos=atof(words[j+3]);
                       } 
                       steerplan->actions[jj].activecv[m].type=1;
                   }else {
                      fprintf(fplog,"|--WRONG STEERPLAN SYNTAX: it must be something like this \n");  
                      fprintf(fplog,"#  this is a comment.. just if needeed \n");  
                      fprintf(fplog,"#  Time  CV n Kspring  Position(wildcard * admitted)  (type CENTRAL/POSITIVE/NEGATIVE, default CENTRAL) CV m Kspring  Position(wildcard * admitted)  (type CENTRAL/POSITIVE/NEGATIVE, default CENTRAL) .... \n");  
                      fprintf(fplog,"#  \n");  
                      fprintf(fplog,"0.0 CV 1 1000.0  *  CV 2 0.0 *  CENTRAL\n");  
                      fprintf(fplog,".....\n");  
                      fprintf(fplog,".....\n");  
                      fprintf(fplog,"\n");  
                      plumed_error("PluMeD dead with errors: check log file"); 
                   }
              };
            }
       jj++; //action counter
    }
    fclose(fp);
    fprintf(fplog,"|--LOADED STEERPLAN: \n");  
    for(i=0;i<steerplan->totstages;i++){
       fprintf(fplog,"|- TIME: %12.6f  ",steerplan->actions[i].t);  
       for(j=0;j<steerplan->ncvs;j++){
         if(steerplan->actions[i].activecv[j].wildcardpos){
          sprintf(string,"*");
         }else {
          sprintf(string,"%12.6f",steerplan->actions[i].activecv[j].pos);
         }
         if(steerplan->actions[i].activecv[j].wildcardk){
          sprintf(string2,"*");
         }else {
          sprintf(string2,"%12.6f",steerplan->actions[i].activecv[j].k);
         }
         if(steerplan->actions[i].activecv[j].type==1){
          sprintf(string3,"CENTRAL ");
         }else if(steerplan->actions[i].activecv[j].type==2){ 
          sprintf(string3,"POSITIVE");
         }else if(steerplan->actions[i].activecv[j].type==3){ 
          sprintf(string3,"NEGATIVE");
         }
         if(steerplan->actions[i].activecv[j].k<0){
              fprintf(fplog," CV %3d SKIPPING THIS STAGE... ",(steerplan->actions[i].activecv[j].ncv+1));
         }else{
              fprintf(fplog," CV %3d KAPPA %12s POS %12s TYPE %8s",steerplan->actions[i].activecv[j].ncv+1,string2,string,string3);
         }
       }
       fprintf(fplog,"\n");  
    }
    // setup colvar header file 
    // 
    // the output should be like: 
    //  STP CV n X val K val T val CV m X val K val T val X val K val     
    sprintf(steerplan->log," XX "); 
    for(i=0;i<steerplan->ncvs;i++){
       int mycv; 
       mycv=steerplan->actions[0].activecv[i].ncv;
       char cvlog[100];
       sprintf(cvlog," XX STPCV%d XX STPX%d XX STPK%d XX STPT%d ",mycv+1,mycv+1,mycv+1,mycv+1);
       strcat(steerplan->log,cvlog);
    }
    
    fprintf(fplog,"|--END LOADING STEERPLAN\n");  
    return 1;
};
