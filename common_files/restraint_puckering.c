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

/*
 * Notation follows
 * M. Sega, E. Autieri and F. Pederiva
 * On the Calculation of Puckering Free Energy Surfaces
 * J. Chem. Phys. 2009 vol. 130 (22) pp. 225102
 *
 */
#define Delta(x,y) ((x)==(y)?5./6.:-1./6.)
#define SQRT3_2 .86602540378443864676

// TRIPLE PRODUCT  a . (bxc)
real PREFIX tprod (rvec a,rvec b, rvec c ) {
	rvec vp;
	oprod(b,c,vp);
	return iprod(a,vp);
}

// WEIGHTS FOR SUMMATION OVER j
// weights for R' and R'' definitions
real puckW[] = {0., SQRT3_2,SQRT3_2,0.,-SQRT3_2,-SQRT3_2};
real puckV[] = {1., 0.5, -0.5 , -1., -0.5, 0.5};
// weights for puckering coordinate definitions
real puckBIGW[] = {0., SQRT3_2,-SQRT3_2,0.,SQRT3_2,-SQRT3_2};
real puckBIGV[] = {1.,-0.5,-0.5,1.,-0.5,-0.5};
real puckSIGN[] = {1.,-1.,1.,-1.,1.,-1.}; // (-1)^(j-1) term in summation

// GENERATE THE VECTORS R_j
// REFERRED TO THE GEOMETRICAL CENTER OF THE RING Rcm
// FROM THE ORIGINAL COORDIANTE VECTORS r_j
real PREFIX generate_R(int i_c,struct mtd_data_s *mtd_data,rvec* R,rvec Rp,rvec Rdp,rvec RpxRdp){
	rvec Rcm,temp;
	real tempr;
	int i,j,alpha,iat,prev;
	rvec folded[6];
	Rcm[0]=Rcm[1]=Rcm[2]=0.0;
	// vector Rcm definition (geometrical center of the ring)
	// position of the first atom of the ring
        iat = colvar.cvatoms[i_c][0];
        folded[0][0]=mtd_data->pos[iat][0];
        folded[0][1]=mtd_data->pos[iat][1];
        folded[0][2]=mtd_data->pos[iat][2];
	// position of the other 5 atom with the minimal image convention
	for(j=1;j<6;j++){ 
            iat = colvar.cvatoms[i_c][j];
            prev = colvar.cvatoms[i_c][j-1];
	    minimal_image(mtd_data->pos[iat],mtd_data->pos[prev],&tempr,temp);
            for(alpha=0;alpha<3;alpha++){
	       folded[j][alpha]=folded[j-1][alpha]+temp[alpha];
	    }
	}
	// definition of the Rcm vector (NOT NORMALIZED!!!)
	for(j=0;j<6;j++){ 
		Rcm[0]+=folded[j][0]; 
		Rcm[1]+=folded[j][1]; 
		Rcm[2]+=folded[j][2]; 
	}
	// vectors R_j definition (ring element's coordinate respect to the center)
	for(j=0;j<6;j++) {
		R[j][0] = folded[j][0]-Rcm[0]/6.; // NOTE: here the Rcm vector is normalized properly
		R[j][1] = folded[j][1]-Rcm[1]/6.; // component by component
		R[j][2] = folded[j][2]-Rcm[2]/6.;
	}  
	// vectors R' and R'' definitions
	Rp[0] = Rp[1] = Rp[2] = 0.;
	Rdp[0] = Rdp[1] = Rdp[2] = 0.;
	for(j=0;j<6;j++){
           for(alpha=0;alpha<3;alpha++){
	      Rp[alpha] += puckW[j] * R[j][alpha];
	      Rdp[alpha] += puckV[j] * R[j][alpha];
	   }
	}
	// vector product R'xR'' definition
	oprod(Rp,Rdp,RpxRdp);
	// returns |R'xR''|
	return  sqrt(iprod(RpxRdp,RpxRdp)); 
}

// PROJECTION z_j OF THE R_j VECTOR
real PREFIX puckering_Zeta (rvec R,rvec RpxRdp,real mod) {
	real tmp;
	tmp = iprod(R,RpxRdp);
	return tmp/mod;
}

// GRADIENT nabla_i OF THE PROJECTION z_j
// returns the vector nabla_i z_j
void PREFIX puckering_gradZeta(rvec gradZ,int index_i,int index_j,real* z,rvec* R,rvec Rp, rvec Rdp, rvec RpxRdp,real mod){
	int l,alpha;
	rvec sum1,sum2,RdpRj,RjRp;
	real Ei,Fi;
	real RpRp,RpRdp,RdpRdp,tmp1,tmp2;
	for(alpha=0;alpha<3;alpha++){
		sum1[alpha]=0.;
		sum2[alpha]=0.;
	}
	// coefficients 
	Ei=0.;	Fi=0.;
	for(l=0;l<6;l++){
		Ei += Delta(index_i,l)*puckW[l] ;
		Fi += Delta(index_i,l)*puckV[l] ;
	}
	// scalar and vector products
	RpRp=iprod(Rp,Rp);
	RdpRdp=iprod(Rdp,Rdp);
	RpRdp=iprod(Rp,Rdp);
	oprod(Rdp,R[index_j],RdpRj);
	oprod(R[index_j],Rp,RjRp);
	// building blocks
	for(alpha=0;alpha<3;alpha++){
		sum1[alpha] = 2.*Ei*( Rp[alpha]*RdpRdp - Rdp[alpha]*RpRdp ) +2.*Fi*( Rdp[alpha]*RpRp - Rp[alpha]*RpRdp);
		sum2[alpha] = Delta(index_i,index_j)*RpxRdp[alpha] + Ei*RdpRj[alpha] + Fi*RjRp[alpha];
	}
	// effective computation of grad_i z_j
        for(alpha=0;alpha<3;alpha++)
	   gradZ[alpha] =-z[index_j]*sum1[alpha]/(2.*mod*mod) +  sum2[alpha]/mod; 
}

// PUCKERING COORDINATE Q (total puckering amplitude)
real PREFIX puckering_Q(real *z){
	real QQ=0.;
	int i;
	for(i=0;i<6;i++) QQ += z[i]*z[i];
	return sqrt(QQ);
}

// GRADIENT OF PUCKERING COORDINATE Q
void PREFIX puckering_gradQ(rvec gradQ, real* z, real Q,rvec*R, rvec Rp, rvec Rdp, rvec RpxRdp,real mod,int index){
	int alpha,j;
	rvec gQ,gzj;
	gQ[0]=gQ[1]=gQ[2]=0.;

	for(j=0;j<6;j++){
		puckering_gradZeta(gzj,index,j,z,R,Rp,Rdp,RpxRdp,mod); 
		for(alpha=0;alpha<3;alpha++){
			gQ[alpha] += z[j] * gzj[alpha];
		}
	}
	for(alpha=0;alpha<3;alpha++){
		gQ[alpha] /= Q;
		gradQ[alpha] = gQ[alpha];
	}
}

// PUCKERING COORDINATE phi
real PREFIX puckering_phi(real *z){
	real tempphi,phi, A, B;
	int j;
	A=0.;	B=0.;
	for(j=0;j<6;j++){
		A += puckBIGW[j]*z[j];
		B += puckBIGV[j]*z[j];
	}
	// disambiguation of the phi angle
	// and its definition in the (CICLIC) range [0,2pi)
	tempphi=atan(-A/B);
	if(B >= 0.){
		if(tempphi >= 0.){
			phi = tempphi;
		}else{	phi = tempphi + 2.*M_PI;}
	}else{ phi = tempphi + M_PI; }
	return phi;
}

// GRADIENT OF PUCKERING COORDINATE phi
void PREFIX puckering_gradphi(rvec gradphi, real* z,rvec*R, rvec Rp, rvec Rdp, rvec RpxRdp,real mod,int index){
	int j,alpha;
	real A,B;
	rvec gphi, gzj,gA,gB;
	gphi[0]=gphi[1]=gphi[2]=0.;
	A=0.;	B=0.;
	gA[0]=gA[1]=gA[2]=0.;
	gB[0]=gB[1]=gB[2]=0.;
	for(j=0;j<6;j++){
		puckering_gradZeta(gzj,index,j,z,R,Rp,Rdp,RpxRdp,mod);
		A += z[j]*puckBIGW[j];
		B += z[j]*puckBIGV[j];
		for(alpha=0;alpha<3;alpha++){
			gA[alpha] += gzj[alpha]*puckBIGW[j];
			gB[alpha] += gzj[alpha]*puckBIGV[j];
		}
	}
	for(alpha=0;alpha<3;alpha++){
		gphi[alpha] = (A*gB[alpha]-B*gA[alpha])/(A*A+B*B);
		gradphi[alpha]=gphi[alpha];
	}	
}

// PUCKERING COORDINATE theta
real PREFIX puckering_theta(real * z){
	real A,B,C;
	int j;
	A=0.;	B=0.;	C=0.;
	for(j=0;j<6;j++){
		A += puckBIGW[j]*z[j];
		B += puckBIGV[j]*z[j];
		C += puckSIGN[j]*z[j];
	}
	return  acos( C / sqrt(2.*(A*A+B*B) +C*C ) );
}

// GRADIENT OF PUCKERING COORDINATE theta
void PREFIX puckering_gradtheta(rvec gradtheta, real* z,rvec*R,rvec Rp,rvec Rdp,rvec RpxRdp,real mod,int index){
	real A,B,C,Q;
	rvec gtheta, gzj,gA,gB,gC;
	int j,alpha;
	gtheta[0]=gtheta[1]=gtheta[2]=0.;
	A=0.;   B=0.;	C=0.;
	gA[0]=gA[1]=gA[2]=0.;
	gB[0]=gB[1]=gB[2]=0.;
	gC[0]=gC[1]=gC[2]=0.;
	Q=puckering_Q(z);
	for(j=0;j<6;j++){
		puckering_gradZeta(gzj,index,j,z,R,Rp,Rdp,RpxRdp,mod);
		A += z[j]*puckBIGW[j];
		B += z[j]*puckBIGV[j];
		C += z[j]*(j%2?-1:1);
		for(alpha=0;alpha<3;alpha++){
			gA[alpha] += gzj[alpha]*puckBIGW[j];
			gB[alpha] += gzj[alpha]*puckBIGV[j];
			gC[alpha] += gzj[alpha]*(j%2?-1:1);
		}
	}
	for(alpha=0;alpha<3;alpha++){
	  gtheta[alpha] = ( C*( B*gB[alpha] + A*gA[alpha] )/sqrt( A*A+B*B ) - sqrt( A*A+B*B)*gC[alpha] );
	  gtheta[alpha] = gtheta[alpha]/( 3.*sqrt(2.)*Q*Q );
	  gradtheta[alpha]=gtheta[alpha];
	}
}

// FOR METADYNAMICS: called by restaint.c 
// this function selects what puckering coordinate will be used
// due to the contents of the META_INP input file
void PREFIX puckering_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i,j, iat;
  real mod,Q;
  rvec rij, sum1, sum2, RpxRdp,Rp,Rdp,gradQ,gradphi,gradtheta;
  static rvec *R=NULL;
  
  real z[6];
  if(R==NULL) R=(rvec*) malloc(sizeof(rvec)*6);

  // generation of the geometrical center coordinate set R_j
  // and accessories vectors (see above)
  mod = generate_R(i_c,mtd_data,R,Rp,Rdp,RpxRdp);

  // generation of the projection z_j
  for(j=0;j<6;j++){
     z[j]=puckering_Zeta (R[j],RpxRdp,mod);
  }
  // calculation of puckering amplitude (needed for all the puckering coords)
  Q = puckering_Q(z);
  // selection of the puckering coordinate
  switch(colvar.type[i_c]) {
	case 1: // puckering coordinate Q
	     { 
                colvar.ss0[i_c] = Q;
		for(i=0;i<6;i++){
		  puckering_gradQ(gradQ, z, Q, R, Rp, Rdp, RpxRdp, mod, i); 
                  colvar.myder[i_c][i][0] = gradQ[0]; 
                  colvar.myder[i_c][i][1] = gradQ[1]; 
                  colvar.myder[i_c][i][2] = gradQ[2]; 
		}
	     }
	       break;
       case 2: // puckering coordinate phi
	    {
		    real phi = puckering_phi(z);
		    rvec gradphi;
		    colvar.ss0[i_c] = phi;
		    for(i=0;i<6;i++){
			puckering_gradphi(gradphi, z, R, Rp, Rdp, RpxRdp, mod, i);
			colvar.myder[i_c][i][0] = gradphi[0];
			colvar.myder[i_c][i][1] = gradphi[1];
			colvar.myder[i_c][i][2] = gradphi[2];
		    }
	    }
	       break;
       case 3: // puckering coordinate theta
	    {
		    real theta = puckering_theta(z);
		    rvec gradtheta;
		    colvar.ss0[i_c] = theta;
		    for(i=0;i<6;i++){
			puckering_gradtheta(gradtheta,z, R, Rp, Rdp, RpxRdp,mod,i);
			colvar.myder[i_c][i][0] = gradtheta[0];
			colvar.myder[i_c][i][1] = gradtheta[1];
			colvar.myder[i_c][i][2] = gradtheta[2];
		    }
		    
	    }
	       break;
       default: break;
		
   }

}

// ------------------------------------------------------------------------------------------------

int PREFIX read_puckering(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw, iat, j;
  double delta = 0.0;
  char string[400];
  int help;
  help=0;

  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR PUCKERING\n"); help=1;}

  if(colvar.natoms[count]!=6) {fprintf(fplog, "PUCKERING is not able to run with more or less than 6 atom in the ring\n"); help=1;}

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  iw=seek_word(word,"TYPE");
  if(iw>0){
    if(seek_word2(word,"Q",iw)>iw)colvar.type[count] = 1;
    if(seek_word2(word,"PHI",iw)>iw)colvar.type[count] = 2;
    if(seek_word2(word,"THETA",iw)>iw)colvar.type[count] = 3;
  } else{ fprintf(fplog,"|- NEEDED TYPE KEYWORD FOR PUCKERING\n"); help=1;}

  if(help){
          printf("%i \n",help);
          fprintf(fplog, "\n-PUCKERING CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "PUCKERING LIST <g1> SIGMA 0.1 TYPE Q/PHI/THETA\n");
          fprintf(fplog, "         g1->            \n");
          fprintf(fplog, "         six atoms       \n");
          fprintf(fplog, "         g1<-            \n");
          fprintf(fplog, "                         \n");
          plumed_error("PluMeD dead with errors: check log file");
  }

  colvar.type_s[count] = 34;

  if(colvar.type[count]==1){
    fprintf(fplog, "\n%i-PUCKERINGQ; ATOMS INVOLVED: %i; ", count+1, colvar.natoms[count]);
  }else if(colvar.type[count]==2){
    fprintf(fplog, "\n%i-PUCKERINGPHI; ATOMS INVOLVED: %i; ", count+1, colvar.natoms[count]);
  }else if(colvar.type[count]==3){
    fprintf(fplog, "\n%i-PUCKERINGTHETA; ATOMS INVOLVED: %i; ", count+1, colvar.natoms[count]);
  }

  snew(colvar.myder[count], colvar.natoms[count]);

 fprintf(fplog,"|- SET MEMBERS: ");
  for(i=0;i<colvar.natoms[count];i++){fprintf(fplog," %d ",colvar.cvatoms[count][i]+1);if((i+1)%20==0)fprintf(fplog,"\n               ");}fprintf(fplog,"\n\n");

  return colvar.natoms[count];
}
