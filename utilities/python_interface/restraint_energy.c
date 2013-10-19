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

void PREFIX energy_restraint(int i_c, struct mtd_data_s *mtd_data)
{
 
  colvar.ss0[i_c] = mtd_data->energy;

}

//-----------------------------------------------------------------------------------------------------------------

int PREFIX read_energy(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iat, iw;
  double delta = 0.0;

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  fprintf(fplog, "\n%1i-ENERGY ", count+1); 
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 

  colvar.natoms[count]   = mtd_data.natoms;

  snew(colvar.myder[count], colvar.natoms[count]);
  snew(colvar.cvatoms[count], colvar.natoms[count]);

  for(iat=0;iat<colvar.natoms[count];iat++) colvar.cvatoms[count][iat] = iat;   

  return colvar.natoms[count];
}
