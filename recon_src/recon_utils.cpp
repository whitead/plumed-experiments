#ifdef RECONMETAD 
#include "recon_types.h"
#define _GAT_EXTERN
#include "recon_utils.h"

namespace std {
rp::mpiostream pout(std::cout), perr(std::cerr);
};

#if defined (PLUMED_GROMACS3) || defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)

#include "gmx_lapack.h"

//extern "C" {

//void F77_FUNC(dgesvd,DGESVD)(const char *jobu, const char *jobvt, int *m, int *n, double *da, int *lda, double* S, double* U, int* ldu, double* VT, int* ldvt, double *work, int *lwork, int *info);

//}

#else

#ifndef DRIVER
extern "C" {

//int F77_FUNC(dgesvd,DGESVD)(const char* jobu, const char* jobvt, int* m, int* n, double* da, int* lda, double* S, double* U, int* ldu, double* VT, int* ldvt, double *work, int *lwork, int *info);
int F77_FUNC(dsyev,DSYEV)( const char* jobz, const char* uplo, int* n, double* da, int* lda, double* evals, double* work, int* lwork, int* info );

}  
#endif

#endif

namespace rp {

// Calculate the distance between two vectors (include periodicity)
real calcDistance( const rp::rVector& period, const rp::rVector& v1, const rp::rVector& v2 ){
#ifdef DEBUG
   if( v1.size()!=v2.size() ) ERROR("SIZE MISMATCH");
   if( period.size()!=v1.size() ) ERROR("PERIOD SIZE MISMATCH");
#endif
   // Establish difference between vectors
   rVector diff( v1 ); diff-=v2; pbc( diff, period );

  return sqrt( dotProduct(diff, diff ) );
}

// Calculate the distance between two vectors (neglect periodicity)
real calcDistance( const rp::rVector& v1, const rp::rVector& v2 ){
#ifdef DEBUG
   if( v1.size()!=v2.size() ) ERROR("SIZE MISMATCH");
#endif
   // Establish difference between vectors
   rVector diff( v1 ); diff-=v2; 

  return sqrt( dotProduct(diff, diff ) );
}

void pbc( rVector& A, const rVector& period ){
#ifdef DEBUG
   if(A.sz!=period.sz) ERROR("Vector dimension mismatch");
#endif
   real tmp;
   for(index i=0;i<A.sz;++i){tmp=pbc(A[i],period[i]); A[i]=tmp;}
}

real pbc( const real& dx, const real& period ){
   if(period==0.){ return dx; }
   else if( dx > period/2.0 ){ return dx-period; }
   else if( dx <= -period/2.0 ){ return dx+period;}
   else{ return dx;}
}

void matMatMult(const rMatrix& A, const rMatrix& B, rMatrix& C){
#ifdef DEBUG
  if(A.cl!=B.rw) ERROR("Matrix dimensions mismatch");
#endif
  if( A.rw!=C.rw || B.cl!=C.cl ){C.resize( A.rw, B.cl );} C=0.;
  for (index i=0;i<A.rw;++i) for (index j=0;j<B.cl;++j) for (index k=0; k<A.cl; ++k) C(i,j)+=A(i,k)*B(k,j);
}

// Multiply the transpose of A times B 
void tmatMatMult(const rMatrix& A, const rMatrix& B, rMatrix& C){
#ifdef DEBUG
  if(A.rw!=B.rw) ERROR("Matrix dimensions mismatch");
#endif
  if( A.cl!=C.rw || B.cl!=C.cl ){C.resize( A.cl, B.cl );} C=0.;
  for (index i=0;i<A.cl;++i) for (index j=0;j<B.cl;++j) for (index k=0; k<A.rw; ++k) C(i,j)+=A(k,i)*B(k,j);
}

void matTMatMult( const rMatrix& A, const rMatrix& B, rMatrix& C){
#ifdef DEBUG
  if(A.cl!=B.cl) ERROR("Matrix dimensions mismatch");
#endif
  if( A.rw!=C.rw || B.rw!=C.cl ){C.resize( A.rw, B.rw );} C=0.;
  for (index i=0;i<A.rw;++i) for (index j=0;j<B.rw;++j) for (index k=0; k<A.cl; ++k) C(i,j)+=A(i,k)*B(j,k);
}

// Multiply a diagonal matrix by a matrix
void dmatMatMult( const rVector& A, const rMatrix& B, rMatrix& C){
#ifdef DEBUG
  if(A.sz!=B.rw) ERROR("Matrix dimensions mismatch");
#endif

  if( A.sz!=C.rw || B.cl!=C.cl ){C.resize( A.sz, B.cl );} C=0.;
  for (index i=0;i<A.sz;++i) for (index j=0;j<B.cl;++j) C(i,j)=A[i]*B(i,j); 
}

void matDMatMult( const rMatrix& A, const rVector& B, rMatrix& C){
#ifdef DEBUG
  if(A.cl!=B.sz) ERROR("Matrix dimensions mismatch");
#endif
  if( A.rw!=C.rw || B.sz!=C.cl ){ C.resize( A.rw, B.sz ); }
  for (index i=0;i<A.rw;++i) for(index j=0;j<B.sz;++j) C(i,j)=A(i,j)*B[j];
}

real dotProduct(const rVector& A, const rVector& B){
#ifdef DEBUG
  if(A.sz!=B.sz) ERROR("Vector dimensions mismatch");
#endif
  real dp=0.;
  for(index i=0;i<A.sz;++i) dp+=A.data[i]*B.data[i]; 
  return dp;
}

void vecMatMult( const rVector& A, const rMatrix& B, rVector& C ){
#ifdef DEBUG
  if( A.sz!=B.rw ) ERROR("Matrix and vector not compatible "<<A.sz<<" "<<B.rw);
#endif
  if( C.sz!=B.cl ){ C.resize(B.cl); } C=0;
  for(index i=0;i<B.cl;++i) for(index k=0;k<A.sz;++k) C.data[i]+=A.data[k]*B(k,i);
}

void matVecMult( const rMatrix& A, const rVector& B, rVector& C){
#ifdef DEBUG
  if( A.cl!=B.sz ) ERROR("Matrix and vector not compatible "<<A.cl<<" "<<B.sz);
#endif
  if( C.sz!=A.rw ){C.resize(A.rw);} C=0.; 
  for(index i=0;i<A.rw;++i) for(index k=0;k<A.cl;++k) C.data[i]+=A(i,k)*B.data[k] ;
}

// Multiply a vector by the transpose of the input matrix 
void TmatVecMult( const rMatrix& A, const rVector& B, rVector& C){
#ifdef DEBUG
  if( A.rw!=B.sz ) ERROR("Matrix and vector not compatible "<<A.rw<<" "<<B.sz);
#endif

  if( C.sz!=A.cl ){C.resize(A.cl);} C=0.;
  for(index i=0;i<A.cl;++i) for(index k=0;k<A.rw;++k) C.data[i]+=A(k,i)*B.data[k] ;
}

void transpose( const rMatrix& A, rMatrix& B){
  B.resize( A.cl, A.rw );
  for(index i=0;i<A.cl;++i) for(index j=0;j<A.rw;++j) B(i,j)=A(j,i);
}

// Don't compile any lapack stuff if we are compiling with driver
#ifndef DRIVER
//void pseudoInverse( const rMatrix& A, rMatrix& pseudoinverse ){
//
//  double *da=new double[A.size()]; index k=0; 
//  // Transfer the matrix to the local array
//  for (index i=0; i<A.cl; ++i) for (index j=0; j<A.rw; ++j) da[k++]=A(j,i);
//
//  // Create column and row information on the matrix
//  int nsv, nrows, ncols, info; nrows=A.rw; ncols=A.cl; 
//  if(nrows>ncols){nsv=ncols;}else{nsv=nrows;}
//
//  // Create some containers for stuff from single value decomposition
//  double *S=new double[nsv]; double *U=new double[nrows*nrows]; 
//  double *VT=new double[ncols*ncols];
//
//  // This optimizes the size of the work array used in lapack singular value decomposition
//  int lwork=-1; double* work=new double[1];
//  F77_FUNC(dgesvd,DGESVD)( "A", "A", &nrows, &ncols, da, &nrows, S, U, &nrows, VT, &ncols, work, &lwork, &info ); 
//  if(info!=0) ERROR("Return "<<info<<" in work optimization call to dgesvd");
//
//  // Retrieve correct sizes for work and rellocate
//  lwork=(int) work[0]; delete [] work; work=new double[lwork];
// 
//  // This does the singular value decomposition
//  F77_FUNC(dgesvd,DGESVD)( "A", "A", &nrows, &ncols, da, &nrows, S, U, &nrows, VT, &ncols, work, &lwork, &info );
//  if(info!=0) ERROR("Return "<<info<<" in call to dgesvd");
//
//  // Compute the tolerance on the singular values ( machine epsilon * nsv * maximum singular value )
//  real tol; tol=S[0]; for(rp::index i=1;i<nsv;++i){ if( S[i]>tol ){ tol=S[i]; } } tol*=nsv*R_EPS;
//
//  // Get the inverses of the singlular values
//  rMatrix Si( ncols , nrows ); Si=0.0;
//  for(rp::index i=0;i<nsv;++i){ if( S[i]>tol ){ Si(i,i)=1./S[i]; }else{ Si(i,i)=0.0; } }
//
//  // Now extract the other matrices we use to compute the pseudoinverse 
//  // ( N.B. these are the transposes of the matrices output by dgesvd )
//  rMatrix V( ncols, ncols ), UT( nrows, nrows ); 
//  k=0; for(rp::index i=0;i<ncols;++i){ for(index j=0;j<ncols;++j){ V(i,j)=VT[k++]; } }
//  
//  k=0; for(rp::index i=0;i<nrows;++i){ for(index j=0;j<nrows;++j){ UT(i,j)=U[k++]; } } 
//
//  // And now compute the psedoinverse
//  rMatrix tmp( ncols, nrows ); pseudoinverse.resize( ncols, nrows );
//  matMatMult( V, Si, tmp ); matMatMult( tmp, UT, pseudoinverse );
//
//}

void invertMat( const rMatrix& A, rMatrix& inverse){

#ifdef DEBUG
   if( A.cl!=A.rw ) ERROR("trying to invert a non square matrix");
#endif

   //compute inverse by diagonalization, since on some machines dgetrf is weirdly broken
   rMatrix evec, tevec; rVector eval; 
   tevec.resize(A.cl,A.rw);
   diagMat(A,eval,evec);
   for (index i=0; i<A.cl; ++i) eval[i]=1.0/eval[i];

   for (index i=0; i<A.rw; ++i) for (index j=0; j<A.cl; ++j) tevec(i,j)=evec(j,i)*eval[j];
   matMatMult(tevec,evec,inverse);

}

void diagMat( const rMatrix& A, rVector& eigenvals, rMatrix& eigenvecs ){

#ifdef DEBUG
   if( A.cl!=A.rw ) ERROR("routine only works for square matrices");
#endif

   double *da=new double[A.size()]; index k=0; double *evals=new double[A.cl]; 
   // Transfer the matrix to the local array
   for (index i=0; i<A.cl; ++i) for (index j=0; j<A.cl; ++j) da[k++]=A(j,i);
   
   int n=A.cl; int lwork=-1, liwork=-1, m, info, one=1; 
   double *work=new double[A.cl]; int *iwork=new int[A.cl];
   double vl, vu, abstol=0.0; 
   int* isup=new int[2*A.cl]; double *evecs=new double[A.sz];

#if defined (PLUMED_GROMACS3) || defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
   F77_FUNC(dsyevr,DSYEVR)("V", "I", "U", &n, da, &n, &vl, &vu, &one, &n ,
                            &abstol, &m, evals, evecs, &n,
                            isup, work, &lwork, iwork, &liwork, &info);
   // Retrieve correct size for liwork and reallocate (this array is not needed when we use dsyev)
   liwork=iwork[0]; delete [] iwork; iwork=new int[liwork];
#else
   F77_FUNC(dsyev,DSYEV)("V", "U", &n, da, &n, evals, work, &lwork, &info );
#endif
   if (info!=0) ERROR("Return "<<info<<" in work optimization call to dsyevd.");

   // Retrieve correct sizes for work and iwork then reallocate
   lwork=(int) work[0]; delete [] work; work=new double[lwork];

#if defined (PLUMED_GROMACS3) || defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
   F77_FUNC(dsyevr,DSYEVR)("V", "I", "U", &n, da, &n, &vl, &vu, &one, &n ,
                            &abstol, &m, evals, evecs, &n,
                            isup, work, &lwork, iwork, &liwork, &info);
#else
   F77_FUNC(dsyev,DSYEV)("V", "U", &n, da, &n, evals, work, &lwork, &info );
#endif
   if (info!=0) ERROR("Return "<<info<<" in call to dsyevd.");

   // Transfer the eigenvalues and eigenvectors to the output 
   eigenvals.resize( A.cl ); eigenvecs.resize( A.cl, A.rw ); k=0;
   for(index i=0;i<A.cl;++i){
      eigenvals[i]=evals[i];
      // N.B. For ease of producing projectors we store the eigenvectors
      // ROW-WISE in the eigenvectors matrix.  The first index is the 
      // eigenvector number and the second the component
#if defined (PLUMED_GROMACS3) || defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
      for(index j=0;j<A.rw;++j){ eigenvecs(i,j)=evecs[k++]; }
#else
      for(index j=0;j<A.rw;++j){ eigenvecs(i,j)=da[k++]; }
#endif
   }
 
   // Deallocate all the memory used by the various arrays
   delete[] da; delete [] work; delete [] evals; delete[] evecs; delete [] iwork; delete [] isup;
   
   return;
}

real logdet( const rMatrix& A ){
#ifdef DEBUG
   if( A.cl!=A.rw ) ERROR("routine only works for square matrices");
#endif

   double *da=new double[A.sz]; index k=0; double *evals=new double[A.cl];
   // Transfer the matrix to the local array
   for (index i=0; i<A.cl; ++i) for (index j=0; j<A.rw; ++j) da[k++]=A(j,i);
 
   int n=A.cl; int lwork=-1, liwork=-1, info, m, one=1;
   double *work=new double[A.cl]; int *iwork=new int[A.cl];
   double vl, vu, abstol=0.0;
   int* isup=new int[2*A.cl]; double *evecs=new double[A.sz];
#if defined (PLUMED_GROMACS3) || defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
   F77_FUNC(dsyevr,DSYEVR)("N", "I", "U", &n, da, &n, &vl, &vu, &one, &n ,
                            &abstol, &m, evals, evecs, &n,
                            isup, work, &lwork, iwork, &liwork, &info);
   liwork=iwork[0]; delete [] iwork; iwork=new int[liwork];
#else
   F77_FUNC(dsyev,DSYEV)("N", "U", &n, da, &n, evals, work, &lwork, &info );
#endif
   if (info!=0) ERROR("Return "<<info<<" in work optimization call to dsyevd.");

   // Retrieve correct sizes for work and iwork then reallocate
   lwork=(int) work[0]; delete [] work; work=new double[lwork];

#if defined (PLUMED_GROMACS3) || defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)
   F77_FUNC(dsyevr,DSYEVR)("N", "I", "U", &n, da, &n, &vl, &vu, &one, &n ,
                            &abstol, &m, evals, evecs, &n,
                            isup, work, &lwork, iwork, &liwork, &info);
#else
   F77_FUNC(dsyev,DSYEV)("N", "U", &n, da, &n, evals, work, &lwork, &info );
#endif
   if (info!=0) ERROR("Return "<<info<<" in call to dsyevd.");
   
   // Transfer the eigenvalues and eigenvectors to the output 
   real det=0; for(index i=0;i<A.cl;i++){ det+=log(evals[i]); }  
   
   // Deallocate all the memory used by the various arrays
   delete[] da; delete [] work; delete [] evals; delete[] evecs; delete [] iwork; delete [] isup;

   return det;
}
#endif

// Does a cholesky decomposition of a matrix
void cholesky( const rMatrix& A, rMatrix& B ){

#ifdef DEBUG
   if( A.cl!=A.rw ) ERROR("routine only works for square matrices");
#endif

   rMatrix L(A.rw,A.rw,0.0); rVector D(A.rw,0.0);
   for(index i=0; i<A.rw; ++i){
      L(i,i)=1.;
      for (index j=0; j<i; ++j){
         L(i,j)=A(i,j);
         for (index k=0; k<j; ++k) L(i,j)-=L(i,k)*L(j,k)*D[k];
         if (D[j]!=0.) L(i,j)/=D[j]; else L(i,j)=0.;
      }
      D[i]=A(i,i);
      for (index k=0; k<i; ++k) D[i]-=L(i,k)*L(i,k)*D[k];
   }

   for(index i=0; i<A.rw; ++i) D[i]=(D[i]>0.?sqrt(D[i]):0.);
   B.resize(A.rw,A.rw); 
   for(index i=0; i<A.sz; ++i) B.data[i]=0.;
   for(index i=0; i<A.rw; ++i) for(index j=0; j<=i; ++j) B(i,j)+=L(i,j)*D[j];
}

// Solves L.y = b where L is a lower triangular matrix that is the solution of a Cholsky 
// decomposition
void chol_elsolve( const rMatrix &M, const rVector &b, rVector &y ){

#ifdef DEBUG
   if( M.cl!=M.rw ) ERROR("routine only works for square matrices");
   if( M(0,1)!=0.0 ) ERROR("matrix is not upper triangular")
   if( b.sz!=M.rw ) ERROR("size mismatch");
#endif

   if(y.sz!=b.sz){ y.resize(b.sz); }

   for(index i=0;i<b.sz;++i){
      y[i]=b[i];
      for(index j=0;j<i;++j) y[i]-=M(i,j)*y[j];
      y[i]*=1.0/M(i,i);
   }
}

void randNDimGauss( const rMatrix &M, const rVector &center, rVector &b){

#ifdef DEBUG
   if( M.cl!=M.rw ) ERROR("routine only works for square matrices");
   if( M(0,1)!=0.0 ) ERROR("matrix is not upper triangular did you forget to cholesky decompose?")
   if( center.sz!=M.rw ) ERROR("dimension mismatch");
#endif

   rVector tmp(M.rw);
   if(tmp.sz!=b.sz){b.resize(tmp.sz);}

   // Generate a set of 1D Gaussian variates
   for (index i=0;i<tmp.sz;++i) tmp[i]=randGauss();
   // A point distributed around the origin 
   matVecMult( M,tmp,b ); 
   // Add the center of the distribution
   b+=center;
}

real randGauss(){
  real v1,v2,rsq,fac;
  static index iset=0; static real gset;
  if(iset==0){
     for(index i=0;i<100;++i){
       v1=2.*drand48()-1.;
       v2=2.*drand48()-1.;
       rsq=v1*v1+v2*v2;
       if (rsq<1. && rsq!=0.){ break; }
     }

     fac=sqrt(-2.*log(rsq)/rsq);
     gset=v1*fac; iset=1;
     return v2*fac;
  }
  else{
     iset=0;
  }

  return gset;
}

real introot(const real& f, const index& n){
     if (f>0) return pow(f,1./n);
     return pow(fabs(f),1./n)*(n%2==0?1.:-1.);
}

real intpow(const real& f, const index& n){
     if (f>0) return pow(f,1.*n);
     return pow(fabs(f),1.*n)*(n%2==0?1.:-1.);
}

std::ostream& operator<<(std::ostream& ostr,const rVector& vec){
   for(index i=0;i<vec.sz;++i) ostr<<vec.data[i]<<" ";
   return ostr;
}

std::istream& operator>>(std::istream& istr,rVector& vec){
   for(index i=0;i<vec.sz;++i){ 
      if( ( istr>>vec.data[i] ).fail() ) ERROR("SOMETHING HAS RONE WRONG IN READING VECTOR");
   }
   return istr;
}

std::ostream& operator<<(std::ostream& ostr,const rMatrix& mat){
   for(index i=0;i<mat.sz;++i) ostr<<mat.data[i]<<" ";
   return ostr;
}

std::istream& operator>>(std::istream& istr,rMatrix& mat){
   for(index i=0;i<mat.sz;++i){
      if( ( istr>>mat.data[i] ).fail() ) ERROR("SOMETHING HAS RONE WRONG IN READING MATRIX"); 
   }
   return istr;
}

char* itoa(int num){
    /* log10(num) gives the number of digits; + 1 for the null terminator */
    int size;
    if(num==0){size=3;}
    else if(num<0){size=log10(-float(num))+2;}
    else{size = log10(float(num)) + 2;}
    char *x =(char *)malloc(size*sizeof(char));
    snprintf(x, size, "%d", num);
    //printf("%d %d\n", num,size);
    return x;
};

real rVector::moment( const index& n ) const {
   real mom; mom=0;
   for(index i=0;i<sz;++i) mom+=intpow(data[i],n);
   mom/=double(sz);
   return introot(mom,n);
}

real rVector::central_moment( const index& n ) const {
   real mean; mean=0;
   for(index i=0;i<sz;++i) mean+=data[i];
   mean/=double(sz);
   real mom; mom=0;
   for(index i=0;i<sz;++i) mom+=intpow( (data[i]-mean),n );  
   mom/=double(sz);
   return mom;             //introot(mom,n);
}

}; //ends namespace rp
#endif
