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

void PREFIX pca_restraint(int i_c, struct mtd_data_s *mtd_data)
{
	int i,j,k,l,n,iat,cvnatom;
        struct rmsd_inpack inpack;
        struct rmsd_mini_outpack outpack;
	rvec rrpfit;
	real cm0,cm1,cm2;
	real cv;
	real dr000,dr010,dr020,dr100,dr110,dr120,dr200,dr210,dr220;
	real dr001,dr011,dr021,dr101,dr111,dr121,dr201,dr211,dr221;
	real dr002,dr012,dr022,dr102,dr112,dr122,dr202,dr212,dr222;
	real d00n,d01n,d02n,d10n,d11n,d12n,d20n,d21n,d22n;
	real r1p0,r1p1,r1p2,sc0,sc1,sc2,myinc0,myinc1,myinc2;
	//NB: I already checked that the atom used are less then MAXATOMS_RMSD
	real egvx[MAXATOMS_RMSD],egvy[MAXATOMS_RMSD],egvz[MAXATOMS_RMSD];



	// **********************************
	// Prepare the system: 
	// - find rotation matrix
	// - move CM in origin
	// **********************************
	if (colvar.pca_align[i_c]==1) {
		// I have to align the structure to the reference (i.e., find rotation matrix and apply it)

		// Fills the inpack structure with the coordinates of the reference (r0) and the current frame (r1)
		inpack.natoms = colvar.pca_align_atoms;
		inpack.totmass = colvar.pca_align_totmass;

		if (!(mtd_data->istep%colvar.pca_upstride[i_c])) {

			//here every colvar.pca_upstride[i_c] steps or every step if UPSTRIDE is not used
			for (i=0;i<colvar.pca_align_atoms;i++){
				iat = colvar.pca_align_list[i];

				inpack.mass[i] = colvar.pca_align_mass[i];	// mass-weighted alignment (but all masses are set to 1)
				//NB: the coordinates have already been translated with CM in the origin upon initialization
				inpack.r0[0][i] = colvar.pca_align_coord[i][0];
				inpack.r0[1][i] = colvar.pca_align_coord[i][1];
				inpack.r0[2][i] = colvar.pca_align_coord[i][2];

				// current conformation
				inpack.r1[0][i] = mtd_data->pos[iat][0];	
				inpack.r1[1][i] = mtd_data->pos[iat][1];
				inpack.r1[2][i] = mtd_data->pos[iat][2];

//				if (logical.debug) printf("%d iat=%d mass=%lf\tr0(% lf % lf % lf)\tr1(% lf % lf % lf)\n",i,iat,inpack.mass[i],inpack.r0[0][i],inpack.r0[1][i],inpack.r0[2][i],inpack.r1[0][i],inpack.r1[1][i],inpack.r1[2][i]);
			}

			// moves r0 and r1 in r0p and r1p with CM in the origin
			// actually r0 is only copied in r0p (iopt=6) since it has already its CM in origin
			// gives in d[][] the rotation matrix and in dd_r1[][][][] its derivatives (simple=0)
			// NB: r0p and r1p are multiplied by the atom mass
			rmsd_mini_pack(inpack,&outpack,6,0,1);


			if (colvar.pca_upstride[i_c]>1) {
				// if I use UPSTRIDE > 1, then save matrix and derivatives in old_* 
				for (i=0;i<3;i++) 
					for (j=0;j<3;j++)  {
						colvar.old_d[i][j] = outpack.d[i][j];
						for (k=0;k<3;k++)
							for (n=0;n<colvar.natoms[i_c];n++) colvar.old_dd_dr1[i][j][k][n] = outpack.dd_dr1[i][j][k][n];
					}
			}

			// puts the new rotation matrix in static variables for efficiency!
			d00n = outpack.d[0][0];
			d01n = outpack.d[0][1];
			d02n = outpack.d[0][2];
			d10n = outpack.d[1][0];
			d11n = outpack.d[1][1];
			d12n = outpack.d[1][2];
			d20n = outpack.d[2][0];
			d21n = outpack.d[2][1];
			d22n = outpack.d[2][2];

		} else {
			// if I'm here, I'm using UPSTRIDE and it is NOT time to calculate a new matrix:
			// I recover it from the last saved one

			// reset the CM of the current positions and fills outpack.r0p/r1p
			cm0 = cm1 = cm2 = 0.;
			for (i=0;i<colvar.pca_align_atoms;i++) {
				iat = colvar.pca_align_list[i];
				cm0 += mtd_data->pos[iat][0]*colvar.pca_align_mass[i];
				cm1 += mtd_data->pos[iat][1]*colvar.pca_align_mass[i];
				cm2 += mtd_data->pos[iat][2]*colvar.pca_align_mass[i];
			}
			cm0 = cm0 / (real) colvar.pca_align_totmass;
			cm1 = cm1 / (real) colvar.pca_align_totmass;
			cm2 = cm2 / (real) colvar.pca_align_totmass;
			for (i=0;i<colvar.pca_align_atoms;i++) {
				iat = colvar.pca_align_list[i];

				//NB: the coordinates have already been translated with CM in the origin upon initialization
				outpack.r0p[0][i] = colvar.pca_align_coord[i][0];
				outpack.r0p[1][i] = colvar.pca_align_coord[i][1];
				outpack.r0p[2][i] = colvar.pca_align_coord[i][2];

				// current conformation
				outpack.r1p[0][i] = mtd_data->pos[iat][0] - cm0;
				outpack.r1p[1][i] = mtd_data->pos[iat][1] - cm1;
				outpack.r1p[2][i] = mtd_data->pos[iat][2] - cm2;
			}

			// restore the matrix derivatives from old one
			for (i=0;i<3;i++) 
				for (j=0;j<3;j++)  {
					for (k=0;k<3;k++)
						for (n=0;n<colvar.natoms[i_c];n++) outpack.dd_dr1[i][j][k][n] = colvar.old_dd_dr1[i][j][k][n];
				}

			// restore the old rotation matrix in static variables for efficiency!
			d00n = colvar.old_d[0][0];
			d01n = colvar.old_d[0][1];
			d02n = colvar.old_d[0][2];
			d10n = colvar.old_d[1][0];
			d11n = colvar.old_d[1][1];
			d12n = colvar.old_d[1][2];
			d20n = colvar.old_d[2][0];
			d21n = colvar.old_d[2][1];
			d22n = colvar.old_d[2][2];
		}

	} else {

		// I am using the NOALIGN keyword: don't need to find the rotation matrix
		// I only calculate the current CM and if also the DIFF keyword is used, I fill the outpack.r0p with
		// the reference structure for later use
		cm0 = cm1 = cm2 = 0.;
		for (i=0;i<colvar.pca_align_atoms;i++) {
			iat = colvar.pca_align_list[i];

			cm0 += mtd_data->pos[iat][0]*colvar.pca_align_mass[i];
			cm1 += mtd_data->pos[iat][1]*colvar.pca_align_mass[i];
			cm2 += mtd_data->pos[iat][2]*colvar.pca_align_mass[i];

			if (colvar.pca_diff[i_c]==1) {
				//NB: the coordinates have already been translated with CM in the origin upon initialization
				outpack.r0p[0][i] = colvar.pca_align_coord[i][0];
				outpack.r0p[1][i] = colvar.pca_align_coord[i][1];
				outpack.r0p[2][i] = colvar.pca_align_coord[i][2];
			}

		}
		cm0 = cm0 / (real) colvar.pca_align_totmass;
		cm1 = cm1 / (real) colvar.pca_align_totmass;
		cm2 = cm2 / (real) colvar.pca_align_totmass;
	}

//	if (logical.debug) printf("PCA Rotation matrix : % lf\t% lf\t% lf\n                  % lf\t% lf\t% lf\n                  % lf\t% lf\t% lf\n",d00n,d01n,d02n,d10n,d11n,d12n,d20n,d21n,d22n);


	// NB :  colvar.natoms[i_c] == colvar.pca_align_atoms
	cvnatom = colvar.pca_align_atoms;

	// copy the coefficients in 3 vectors for efficiency
	for(n=0; n<cvnatom; n++) {
		egvx[n] = colvar.pcacomp[i_c][n].coeff[0];
		egvy[n] = colvar.pcacomp[i_c][n].coeff[1];
		egvz[n] = colvar.pcacomp[i_c][n].coeff[2];
	}


	// **********************************
	// Calculates new CV value
	// **********************************
	cv = 0.;
	for(n=0; n<cvnatom; n++){
//		if (logical.debug) printf("outpack %d r0(% lf % lf % lf)\tr1(% lf % lf % lf)\n",n,outpack.r0p[0][n],outpack.r0p[1][n],outpack.r0p[2][n],outpack.r1p[0][n],outpack.r1p[1][n],outpack.r1p[2][n]);

		if (colvar.pca_align[i_c]==1) {
			// apply the rotation to the current position (translated with CM in origin)
			rrpfit[0] = d00n * outpack.r1p[0][n] + d01n * outpack.r1p[1][n] + d02n * outpack.r1p[2][n];
			rrpfit[1] = d10n * outpack.r1p[0][n] + d11n * outpack.r1p[1][n] + d12n * outpack.r1p[2][n];
			rrpfit[2] = d20n * outpack.r1p[0][n] + d21n * outpack.r1p[1][n] + d22n * outpack.r1p[2][n];
		} else {
			// I don't have to align to reference, only use the current conformation (translated with CM in origin)
			iat = colvar.pca_align_list[n];
			rrpfit[0] = mtd_data->pos[iat][0] - cm0;
			rrpfit[1] = mtd_data->pos[iat][1] - cm1;
			rrpfit[2] = mtd_data->pos[iat][2] - cm2;
		}

		if (colvar.pca_diff[i_c]==1) {
			// g_anaeig version :
			// project the differences between current roto-translated (=aligned) structure and reference onto the eigenvector
			// NB: since r0p and r1p are returned multiplied by the atom mass, I divide their contribution by the mass
			cv +=((rrpfit[0] - outpack.r0p[0][n]) * egvx[n] \
			    + (rrpfit[1] - outpack.r0p[1][n]) * egvy[n] \
			    + (rrpfit[2] - outpack.r0p[2][n]) * egvz[n]) / colvar.pca_align_mass[n];
		} else {
			// directly project the current roto-translated (=aligned to reference) structure onto the eigenvector
			cv +=(rrpfit[0] * egvx[n] \
			    + rrpfit[1] * egvy[n] \
			    + rrpfit[2] * egvz[n]) / colvar.pca_align_mass[n];
		}
	}


	//returns the value of the CV
	colvar.ss0[i_c] = cv;

//	if (logical.debug) printf("CV = %lf\n",cv);


	// **********************************
	// Calculates CV derivative
	// **********************************
	if (colvar.pca_align[i_c]==0) {

		// If I don't align (ie rotate) the derivative is the egv coeff minus the contribution for the CM displacement
		for(i=0; i<cvnatom; i++) {
			colvar.myder[i_c][i][0] = egvx[i];
			colvar.myder[i_c][i][1] = egvy[i];
			colvar.myder[i_c][i][2] = egvz[i];

			myinc0 = 0;
			myinc1 = 0;
			myinc2 = 0;
			for(j=0;j<cvnatom;j++) {
				myinc0 += egvx[j];
				myinc1 += egvy[j];
				myinc2 += egvz[j];
			}

			colvar.myder[i_c][i][0] -= myinc0/(real) colvar.pca_align_totmass;
			colvar.myder[i_c][i][1] -= myinc1/(real) colvar.pca_align_totmass;
			colvar.myder[i_c][i][2] -= myinc2/(real) colvar.pca_align_totmass;
		}

	} else {

		// If I align (ie rotate) it is more complicated
		for(i=0; i<cvnatom; i++) {
			colvar.myder[i_c][i][0] = d00n*egvx[i] + d10n*egvy[i] + d20n*egvz[i];
			colvar.myder[i_c][i][1] = d01n*egvx[i] + d11n*egvy[i] + d21n*egvz[i];
			colvar.myder[i_c][i][2] = d02n*egvx[i] + d12n*egvy[i] + d22n*egvz[i];
		}

		d00n /= (real) inpack.totmass;
		d01n /= (real) inpack.totmass;
		d02n /= (real) inpack.totmass;
		d10n /= (real) inpack.totmass;
		d11n /= (real) inpack.totmass;
		d12n /= (real) inpack.totmass;
		d20n /= (real) inpack.totmass;
		d21n /= (real) inpack.totmass;
		d22n /= (real) inpack.totmass;
	
		for(i=0; i<cvnatom; i++) {

			dr000 = outpack.dd_dr1[0][0][0][i];
			dr010 = outpack.dd_dr1[0][1][0][i];
			dr020 = outpack.dd_dr1[0][2][0][i];
			dr100 = outpack.dd_dr1[1][0][0][i];
			dr110 = outpack.dd_dr1[1][1][0][i];
			dr120 = outpack.dd_dr1[1][2][0][i];
			dr200 = outpack.dd_dr1[2][0][0][i];
			dr210 = outpack.dd_dr1[2][1][0][i];
			dr220 = outpack.dd_dr1[2][2][0][i];

			dr001 = outpack.dd_dr1[0][0][1][i];
			dr011 = outpack.dd_dr1[0][1][1][i];
			dr021 = outpack.dd_dr1[0][2][1][i];
			dr101 = outpack.dd_dr1[1][0][1][i];
			dr111 = outpack.dd_dr1[1][1][1][i];
			dr121 = outpack.dd_dr1[1][2][1][i];
			dr201 = outpack.dd_dr1[2][0][1][i];
			dr211 = outpack.dd_dr1[2][1][1][i];
			dr221 = outpack.dd_dr1[2][2][1][i];

			dr002 = outpack.dd_dr1[0][0][2][i];
			dr012 = outpack.dd_dr1[0][1][2][i];
			dr022 = outpack.dd_dr1[0][2][2][i];
			dr102 = outpack.dd_dr1[1][0][2][i];
			dr112 = outpack.dd_dr1[1][1][2][i];
			dr122 = outpack.dd_dr1[1][2][2][i];
			dr202 = outpack.dd_dr1[2][0][2][i];
			dr212 = outpack.dd_dr1[2][1][2][i];
			dr222 = outpack.dd_dr1[2][2][2][i];

			myinc0 = 0;
			myinc1 = 0;
			myinc2 = 0;

			for(j=0;j<cvnatom;j++) {

				r1p0 = outpack.r1p[0][j];		// NB: r1p is actually the position*mass
				r1p1 = outpack.r1p[1][j];
				r1p2 = outpack.r1p[2][j];

				sc0=egvx[j];
				sc1=egvy[j];
				sc2=egvz[j];

				// NB: the '-' part comes from the translation of the CM in the origin
				myinc0 += (dr000*r1p0 + dr010*r1p1 + dr020*r1p2 - d00n) * sc0 \
					+ (dr100*r1p0 + dr110*r1p1 + dr120*r1p2 - d10n) * sc1 \
					+ (dr200*r1p0 + dr210*r1p1 + dr220*r1p2 - d20n) * sc2;

				myinc1 += (dr001*r1p0 + dr011*r1p1 + dr021*r1p2 - d01n) * sc0 \
					+ (dr101*r1p0 + dr111*r1p1 + dr121*r1p2 - d11n) * sc1 \
					+ (dr201*r1p0 + dr211*r1p1 + dr221*r1p2 - d21n) * sc2;

				myinc2 += (dr002*r1p0 + dr012*r1p1 + dr022*r1p2 - d02n) * sc0 \
					+ (dr102*r1p0 + dr112*r1p1 + dr122*r1p2 - d12n) * sc1 \
					+ (dr202*r1p0 + dr212*r1p1 + dr222*r1p2 - d22n) * sc2;


			}

			colvar.myder[i_c][i][0] += myinc0;
			colvar.myder[i_c][i][1] += myinc1;
			colvar.myder[i_c][i][2] += myinc2;

/*			if (logical.debug) {
				printf("der[i=%d] x = %lf (myinc0=%e)\n",i,colvar.myder[i_c][i][0],myinc0);
				printf("der[i=%d] y = %lf (myinc1=%e)\n",i,colvar.myder[i_c][i][1],myinc1);
				printf("der[i=%d] z = %lf (myinc2=%e)\n",i,colvar.myder[i_c][i][2],myinc2);
			}
*/
		}

	}
	return;
}

// ------------------------------------------------------------------------------------------------


int PREFIX read_pca(char **word, int count, t_plumed_input *input, FILE *fplog)
{
	int i,n,iw,a1,a2,natom,align_natom,help;
	double x,y,z,coeff,mass;
	char aux[100];
	FILE *fp;

	help = 0;
	natom = 0;
	n = 0;

	colvar.type_s[count] = 42;	// sets the CV type
	//	colvar.cell_pbc[count]=1;	// default is PBC

	colvar.pca_diff[count] = 0;		//default, don't do difference (to do it, use keyword DIFF)
	colvar.pca_align[count] = 1;		//default, align structure to reference befor projection (to avoid alignment, use keyword NOALIGN)
	colvar.pca_upstride[count] = -1;	//default, update the rotation matrix and derivatives every 1 step

	// ***************************************************** EIGENVEC ******************************
	// ***************************************************** EIGENVEC ******************************
	iw = seek_word(word,"EIGENVEC");
	if(iw>=0){   
		fp = fopen(word[iw+1],"r");
		if (fp==NULL) {
			fprintf(fplog,"\n-PCA CV: error opening EIGENVEC file [%s]\nAborting\n\n",word[iw+1]);
			fprintf(stderr,"\nPLUMED ERROR : PCA CV: error opening EIGENVEC file [%s]\nAborting\n\n",word[iw+1]);
			EXIT();
		}

		//do a first loop just to count the number of coefficients. I need it to allocate the appropriate memory
		while(fgets(aux,100,fp)!=NULL) {
			//ignoring comment lines
			if (strncmp(aux,"#",1)) {
				if(sscanf(aux,"%d %lf %lf %lf",&a1,&x,&y,&z)!=4) {
					fprintf(fplog,"\n-PCA CV: error reading EIGENVEC file [%s]: expecting 4 columns: atomid x y z\nAborting\n\n",word[iw+1]);
					fprintf(stderr,"\nPLUMED ERROR : PCA CV: error reading EIGENVEC file [%s]: expecting 4 columns: atomid x y z\nAborting\n\n",word[iw+1]);
					EXIT();
				}
				natom++;	//number of atoms
			}
		}

		//allocate structure pcacomp
		snew(colvar.pcacomp[count],natom);

		//allocate memory for atoms, coefficients and derivatives
		snew(colvar.cvatoms[count], natom);	//will contain the atom id
		snew(colvar.myder[count], natom);	//will contain the derivatives

		//saves the number of atoms
		colvar.natoms[count] = natom;

		rewind(fp);
		while(fgets(aux,100,fp)!=NULL) {
			//ignoring comment lines
			if (strncmp(aux,"#",1)) {
				if(sscanf(aux,"%d %lf %lf %lf",&a1,&x,&y,&z)!=4) {
					fprintf(fplog,"\n-PCA CV: error reading EIGENVEC file [%s]: expecting 4 columns: atomid x y z\nAborting\n\n",word[iw+1]);
					fprintf(stderr,"\nPLUMED ERROR : PCA CV: error reading EIGENVEC file [%s]: expecting 4 columns: atomid x y z\nAborting\n\n",word[iw+1]);
					EXIT();
				}
				colvar.pcacomp[count][n].a1 = a1-1;
				colvar.pcacomp[count][n].coeff[0] = x;
				colvar.pcacomp[count][n].coeff[1] = y;
				colvar.pcacomp[count][n].coeff[2] = z;
				n++;
			}
		}

		fclose(fp);

		if (logical.debug) { 
			fprintf(fplog,"\n%1i-PCA CV: correctly read %d atoms and coefficients from input EIGENVEC file [%s]",count+1,natom,word[iw+1]);
			fflush(NULL);
		}
	} else { fprintf(fplog,"|- NEEDED EIGENVEC KEYWORD FOR PCA CV\n"); help=1; }


	// ***************************************************** FRAME ******************************
	// ***************************************************** FRAME ******************************
	align_natom = n = 0;
	iw=seek_word(word,"FRAME");
	if(iw>=0){   
		fp = fopen(word[iw+1],"r");
		if (fp==NULL) {
			fprintf(fplog,"\n-PCA CV: error opening FRAME file [%s]\nAborting\n\n",word[iw+1]);
			fprintf(stderr,"\nPLUMED ERROR : PCA CV: error opening FRAME file [%s]\nAborting\n\n",word[iw+1]);
			EXIT();
		}

		//do a first loop just to count the number of atoms. I need it to allocate the appropriate memory
		while(fgets(aux,100,fp)!=NULL) {
			//ignoring comment lines
			if (strncmp(aux,"#",1)) {
				if(sscanf(aux,"%d %lf %lf %lf",&a1,&x,&y,&z)!=4) {
					fprintf(fplog,"\n-PCA CV: error reading FRAME file [%s]: expecting 4 columns: atomid x y z\nAborting",word[iw+1]);
					fprintf(stderr,"\nPLUMED ERROR : PCA CV: error reading FRAME file [%s]: expecting 4 columns: atomid x y z\nAborting",word[iw+1]);
					EXIT();
				}
				// DEVEL : expect also a mass column
				/*                        	if(sscanf(aux,"%d %lf %lf %lf %lf",&a1,&x,&y,&z,&mass)!=5) {
								fprintf(fplog,"\n-PCA CV: error reading FRAME file [%s]: expecting 5 columns: atomid x y z mass\nAborting",word[iw+1]);
								fprintf(stderr,"\nPLUMED ERROR : PCA CV: error reading FRAME file [%s]: expecting 5 columns: atomid x y z mass\nAborting",word[iw+1]);
								EXIT();
								}
				 */
				align_natom++;	//number of atoms
			}
		}

		colvar.pca_align_atoms = align_natom;
		snew(colvar.pca_align_list, align_natom);	// contain the list of atoms of the reference structure used to align
		snew(colvar.pca_align_coord, align_natom);	// contain the xyz coordinates of the reference structure 
		snew(colvar.pca_align_mass, align_natom);	// contain the list of atoms of the reference structure used to align
		colvar.pca_align_totmass = 0.;

		rewind(fp);
		while(fgets(aux,100,fp)!=NULL) {
			//ignoring comment lines
			if (strncmp(aux,"#",1)) {
				if(sscanf(aux,"%d %lf %lf %lf",&a1,&x,&y,&z)!=4) {
					fprintf(fplog,"\n-PCA CV: error reading FRAME file [%s]: expecting 4 columns: atomid x y z\nAborting",word[iw+1]);
					fprintf(stderr,"\nPLUMED ERROR : PCA CV: error reading FRAME file [%s]: expecting 4 columns: atomid x y z\nAborting",word[iw+1]);
					EXIT();
				}

				mass = 1.0;	// forces mass = 1, the use of a mass is under development

				// DEVEL : expect also a mass column
				/*                        	if(sscanf(aux,"%d %lf %lf %lf %lf",&a1,&x,&y,&z,&mass)!=5) {
								fprintf(fplog,"\n-PCA CV: error reading FRAME file [%s]: expecting 5 columns: atomid x y z mass\nAborting",word[iw+1]);
								fprintf(stderr,"\nPLUMED ERROR : PCA CV: error reading FRAME file [%s]: expecting 5 columns: atomid x y z mass\nAborting",word[iw+1]);
								EXIT();
								}
				 */
				if (!(mass>0)) plumed_error("null mass not allowed for atoms in reference frame file");
				colvar.pca_align_list[n] = a1-1;
				colvar.cvatoms[count][n] = a1-1;
				colvar.pca_align_coord[n][0] = x;	//NB: I leave unaltered the units
				colvar.pca_align_coord[n][1] = y;
				colvar.pca_align_coord[n][2] = z;
				colvar.pca_align_mass[n] = mass;
				colvar.pca_align_totmass += mass;
				n++;
			}
		}
		fclose(fp);

		if (colvar.pca_align_totmass<0.000001) {
			fprintf(fplog,"\n-PCA CV: ERROR total mass is : %lf\n-PCA CV: if you do not want to mass-weight the alignment to the reference frame, just set to 1 every atom mass.\n",colvar.pca_align_totmass);
			EXIT();
		}

		// centers the atoms with the CM in the origin (so I don't have to do it every step)
		x = y = z = 0.;
		for (i=0;i<natom;i++) {
			x += colvar.pca_align_coord[i][0] * colvar.pca_align_mass[i];
			y += colvar.pca_align_coord[i][1] * colvar.pca_align_mass[i];
			z += colvar.pca_align_coord[i][2] * colvar.pca_align_mass[i];
		}
		x = x / colvar.pca_align_totmass;
		y = y / colvar.pca_align_totmass;
		z = z / colvar.pca_align_totmass;
		if (logical.debug) printf(" CENTERING ref frame (old CM :  % lf\t% lf\t% lf)\n",x,y,z);
		for (i=0;i<natom;i++) {
			colvar.pca_align_coord[i][0] -= x;
			colvar.pca_align_coord[i][1] -= y;
			colvar.pca_align_coord[i][2] -= z;
			if (logical.debug) printf("%i = % lf\t% lf\t% lf\n",i,colvar.pca_align_coord[i][0],colvar.pca_align_coord[i][1],colvar.pca_align_coord[i][2]);
		}

		if (logical.debug) { 
			fprintf(fplog,"\n%1i-PCA CV: correctly read %d atoms and coordinates from input FRAME file [%s] with total mass %lf",count+1,natom,word[iw+1],colvar.pca_align_totmass);
			fflush(NULL);
		}
	} else { fprintf(fplog,"|- NEEDED FRAME KEYWORD FOR PCA CV\n"); help=1; }

	// ***************************************************** DIFF ****************************
	// ***************************************************** DIFF ****************************
	iw=seek_word(word,"DIFF");
	if(iw>=0) colvar.pca_diff[count] = 1;

	// ***************************************************** NOALIGN ****************************
	// ***************************************************** NOALIGN ****************************
	iw=seek_word(word,"NOALIGN");
	if(iw>=0) colvar.pca_align[count] = 0;


	// ***************************************************** UPSTRIDE **************************
	// ***************************************************** UPSTRIDE **************************
	iw=seek_word(word,"UPSTRIDE");		//UNDER DEVEL!!
	if(iw>=0) {
		sscanf(word[iw+1],"%lf", &coeff);
		colvar.pca_upstride[count] = (int) coeff;
	}

	// ***************************************************** SIGMA ******************************
	// ***************************************************** SIGMA ******************************
	iw=seek_word(word,"SIGMA");
	if(iw>=0) {
		sscanf(word[iw+1],"%lf", &coeff);
		colvar.delta_r[count]  = (real) coeff;
	}
	// else { fprintf(fplog,"|- NEEDED SIGMA KEYWORD FOR PCA\n"); help=1; }


	// ********************** I CHECK THAT THE ATOMS ON WHICH I ALIGN ARE THE SAME SET OF THE ATOMS ON WHICH I CALCULATED THE PCA ******
	// ********************** I CHECK THAT THE ATOMS ON WHICH I ALIGN ARE THE SAME SET OF THE ATOMS ON WHICH I CALCULATED THE PCA ******
	// ********************** I CHECK THAT THE ATOMS ON WHICH I ALIGN ARE THE SAME SET OF THE ATOMS ON WHICH I CALCULATED THE PCA ******
	if (colvar.pca_align_atoms != natom) {fprintf(fplog,"ERROR : THE NUMBER OF ATOMS ON WHICH I ALIGN IS DIFFERENT FROM THE NUMBER OF COMPONENTS\n");EXIT();}
	if (colvar.pca_align_atoms > MAXATOMS_RMSD)  {fprintf(fplog,"ERROR : INCREASE MAXATOMS_RMSD IN metadyn.h FILE TO AT LEAST THE NUMBER OF ATOMS USED : %d\n",colvar.pca_align_atoms);EXIT();}
	for (n=0;n<natom;n++) {
		if (colvar.pca_align_list[n] != colvar.pcacomp[count][n].a1 ) {fprintf(fplog,"ERROR : THE ATOMS ON WHICH I ALIGN IS DIFFERENT FROM THE ATOMS ON WHICH I APPLY THE EIGENVECTOR\n");EXIT();}
	}


	if(help){
		fprintf(fplog, "\nPCA CV: WRONG SYNTAX\n");
		fprintf(fplog, "e.g.:     \n");
		fprintf(fplog, "PCA FRAME reference.dat EIGENVEC eigvector.dat [DIFF] [SIGMA s]\n");
		//		fprintf(fplog, "PCA FRAME reference.dat EIGENVEC eigvector.dat [DIFF] [UPSTRIDE n] [SIGMA s]\n");	//DEVEL
		fprintf(fplog, "where\n");
		fprintf(fplog, " reference.dat is a 4-columns file with the reference structure used for alignment in the format: atomid x y z\n");
		fprintf(fplog, " eigvector.dat is a 4-columns file with the eigenvector in the format: atomid x y z\n\n");
		fprintf(fplog, " if DIFF keyword is used, the difference between the current conformation and the reference frame is projected: CV=<(X-ref),evec>\n\n");
		//		fprintf(fplog, " if UPSTRIDE keyword is used, the rotation matrix and derivatives are updated every n steps instead of every step (for performance) NB UNDER TESTING!!!!! \n\n");//DEVEL
		fprintf(fplog, "NB: the alignment is non mass weighted\n");
		fprintf(fplog, "NB: double precision code and ALIGN_ATOMS directive are recommended\n");
		fprintf(fplog, "NB: x y z are intended in engine units\n");
		fprintf(fplog, "NB: in both files, lines beginning with '#' are ignored, no white lines are allowed\n");
		//		fprintf(fplog, "NB: UPSTRIDE should be smaller than HILLS STRIDE\n");		//DEVEL
		plumed_error("PluMed dead with errors: check log file\n");
	}

	fprintf(fplog, "\n%1i-PCA: %d atoms used for alignment and projection",count+1,align_natom);

	if (colvar.pca_upstride[count]>1) fprintf(fplog,"; UPSTRIDE %d\n",colvar.pca_upstride[count]);
	else colvar.pca_upstride[count] = 1;	//set the default to every step

	if (logical.do_hills) fprintf(fplog,"; SIGMA %f\n",colvar.delta_r[count]);
	else fprintf(fplog,"\n");

	if (colvar.pca_align[count]==1) {
		if (colvar.pca_diff[count]==1) {
			//align and diff
			fprintf(fplog,"|- DIFF KEYWORD PRESENT\n");
			fprintf(fplog,"|- ALIGN CURRENT STRUCTURE X TO REFERENCE R AND PROJECT THE DIFFERENCE WITH R : CV = < (X'-R), EIGENVEC >\n");
		} else {
			//align no diff
			fprintf(fplog,"|- ALIGN CURRENT STRUCTURE X TO REFERENCE R AND PROJECT IT : CV = < X', EIGENVEC >\n");
		}
	} else {
		if (colvar.pca_diff[count]==1) {
			//no align, yes diff
			fprintf(fplog,"|- NOALIGN AND DIFF KEYWORDS PRESENT\n");
			fprintf(fplog,"|- NO ALIGNMENT OF CURRENT STRUCTURE X, PROJECTION OF THE DIFFERENCE WITH REFERENCE R : CV = < (X-R), EIGENVEC >\n");
		} else {
			//no align no diff
			fprintf(fplog,"|- NOALIGN KEYWORD PRESENT\n");
			fprintf(fplog,"|- DIRECT PROJECTION OF THE CURRENT STRUCTURE X : CV = < X, EIGENVEC >\n");
		}
	}

	fprintf(fplog,"|- EIGENVECTOR COMPONENTS\n|-  #  \t: ATOM\t X Y Z\n");
	for(i=0;i<natom;i++){
		fprintf(fplog,"|- %3d\t: %d \t % lf % lf % lf\n",i+1,colvar.pcacomp[count][i].a1+1,colvar.pcacomp[count][i].coeff[0],colvar.pcacomp[count][i].coeff[1],colvar.pcacomp[count][i].coeff[2]);
	}
	fflush(NULL);

	return colvar.natoms[count];
}

