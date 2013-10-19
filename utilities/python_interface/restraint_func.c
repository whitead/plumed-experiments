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

void PREFIX func_restraint(int i_c, struct mtd_data_s *mtd_data)
{

int i,ix;
int jcv;
int iatG, iatT;
real ff;

#ifdef HAVE_MATHEVAL
   struct mathfunction_s *myfunc;
   myfunc=&mathfunction[i_c];
   // reset all
   for(i=0;i<colvar.natoms[i_c];i++){ 
   for(ix=0;ix<3;ix++){ 
      colvar.myder[i_c][i][ix]=0; } }

   i=0;
   for (i = 0; i < myfunc->count; i++){
       // collect all the needed vals
       myfunc->valcvs[i]=colvar.ss0[ myfunc->indcvs[i] ]; 
       // calculate the function
   }
   colvar.ss0[i_c]=evaluator_evaluate(myfunc->f,myfunc->count, myfunc->names,myfunc->valcvs); 
   // now calculate the derivatives
   iatG=0;
   for (jcv = 0; jcv < myfunc->count; jcv++){
	  ff=evaluator_evaluate(myfunc->f_prim[jcv],myfunc->count, myfunc->names,myfunc->valcvs); 
      for(iatT=0;iatT<colvar.natoms[myfunc->indcvs[jcv]];iatT++){
         for(ix=0;ix<3;ix++) colvar.myder[i_c][iatG][ix] += ff*colvar.myder[myfunc->indcvs[jcv]][iatT][ix];
         iatG++;
      }
    }

#else
   fprintf(mtd_data->fplog, "FUNCTION OF CVS NOT SUPPORTED: COMPILE WITH LIBMATHEVAL \n");
   plumed_error("PluMeD dead with errors: check log file");
   EXIT(); 
#endif

return;
}



int PREFIX read_func(char **word, int count, t_plumed_input *input,int *iline, int *nw ,FILE *fplog)
{
char *tok,tok2;
int iw, i, j, k;
double delta=0.0;
int nterms;
int help;
int go;
int done;
help=0;
done=0;
int jcv;
char temp[300];
struct mathfunction_s *myfunc;
#ifdef HAVE_MATHEVAL
myfunc=&mathfunction[count];

   go=0; 
   //fprintf(fplog,"NW IS %d \n",*nw);
   // allow spacing: from first quote up to last quote
   for(i=0;i<(*nw);i++){
       // find the quotes 
       tok=strstr(word[i],"\"");
       if(tok!=NULL){
           if(go==1){ go=0; 
                 // fprintf(fplog,"LASTTOK IS %s WORD IS %s\n",tok,word[i]);
                  strcat(myfunc->fline,word[i]);
                  strcat(myfunc->fline,"\0");done=1;}
           else     { go=1; 
                 // fprintf(fplog,"FIRSTTOK IS %s WORD IS %s\n",tok,word[i]); 
                  strcpy(myfunc->fline,word[i]) ; };
       };
       if(go==1 && tok==NULL ){  
                 //fprintf(fplog,"WORD IS %s\n",word[i])  ; 
                   strcat(myfunc->fline,word[i])  ;}
   } 
   if(done==0)help=1;
   if(help){
    fprintf(fplog, "\n|-FUNCTION CV: WRONG SYNTAX\n");
    fprintf(fplog, "e.g.:     \n");
    fprintf(fplog, "      FUNCTION \" CV_1+3*CV_2 \"  {SIGMA 1.0}\n");
    plumed_error("PluMeD dead with errors: check log file");
   }


   //fprintf(fplog,"THE FUNCTION IS %s \n",myfunc->fline); 
   // now substitute
   i=0;
   strcpy(temp,"");
   while(myfunc->fline[i] !='\0'){
     if( strncmp(&myfunc->fline[i],"\"",1)!=0){ strncat(temp,&myfunc->fline[i],1);}
     i++;
   }
   //EXIT();
   strcat(temp,"\0"); 
   strcpy(myfunc->fline,temp);
   fprintf(fplog, "\n%1i-FUNCTION: %s \n", count+1,myfunc->fline);
 
   myfunc->f = evaluator_create (myfunc->fline);
   evaluator_get_variables (myfunc->f, &myfunc->names, &myfunc->count);
   // allocate the primes 300 chars should be enough
   myfunc->f_prim=(void **)malloc(myfunc->count*sizeof(void *));
   myfunc->indcvs=(int *)malloc(myfunc->count*sizeof(int));
   myfunc->valcvs=(real *)malloc(myfunc->count*sizeof(real));
   fprintf(fplog,"|-VARIABLES IDENTIFIED: \n");
   // parse the names
   k=0;
   colvar.natoms[count]=0;
   for (i = 0; i < myfunc->count; i++){
          fprintf(fplog, "|-CV %d  :   %s \n",i+1, myfunc->names[i]);
          tok=strtok(myfunc->names[i],"CV_");
          myfunc->indcvs[i]=atoi(tok)-1; 
          logical.always[  myfunc->indcvs[i]  ]=1; 
          if(myfunc->indcvs[i]>=count){
              fprintf(fplog,"\n\nCV %d DOES NOT APPEAR TO HAVE BEEN DEFINED YET.\n",jcv+1);
              fprintf(fplog,"ONLY PREVIOUSLY DEFINED CVs CAN BE COMBINED BY FUNC!\n");
              plumed_error("PluMeD dead with errors: check log file");
              EXIT(); // just to be sure!
          };
          //fprintf(fplog, "INDEXCV  :   %d \n",myfunc->indcvs[i]);
          jcv=myfunc->indcvs[i];
          colvar.natoms[count]+=colvar.natoms[jcv]; 
         // fprintf(fplog,"jcv=%d, natoms=%d\n", jcv, colvar.natoms[jcv]);
          srenew(colvar.cvatoms[count],colvar.natoms[count]);
      // Copy atom indexes from jcv into count 
          for(j=0;j<colvar.natoms[jcv];j++){
             colvar.cvatoms[count][k] = colvar.cvatoms[jcv][j]; 
           k++; }
     }

     for (i = 0; i < myfunc->count; i++){
       fprintf (fplog,"|- FIRST DERIVATIVE RESPECT TO %s :  ", myfunc->names[i]);
       myfunc->f_prim[i]=evaluator_derivative(myfunc->f,myfunc->names[i]);
       fprintf (fplog," %s\n", evaluator_get_string (myfunc->f_prim[i]));
     } 

fprintf(fplog,"|%d ATOMS IN THIS CV:\n", colvar.natoms[count]);
for(j=0; j<colvar.natoms[count]; j++){
  fprintf(fplog, "|-(ifunc=%d), (atid=%d)\n",j,colvar.cvatoms[count][j]);
} 
 
colvar.type_s[count]   = 51;
snew(colvar.myder[count], colvar.natoms[count]);

iw=seek_word(word,"SIGMA");
if(iw>=0){ 
  sscanf(word[iw+1],"%lf", &delta);
  colvar.delta_r[count]  = (real) delta; 
}

 if (logical.do_hills) fprintf(fplog,"|- SIGMA %f\n",colvar.delta_r[count]);
 else fprintf(fplog,"\n");

 return colvar.natoms[count];
// EXIT(); 
#else
  fprintf(fplog, "\n-FUNCTION CV: YOU CANNOT DO THAT. COMPILE WITH LIBMATHEVAL \n");
  plumed_error("PluMeD dead with errors: check log file");
  EXIT();
#endif
}


