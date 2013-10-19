#ifndef __rp_utilsh
#define __rp_utilsh

#include "recon_types.h"
#include <valarray>
#include <stdlib.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>

namespace rp {
class rMatrix;
class rVector;

// Calculates the nth root of a number
real introot( const real& f, const index& n);
// Calculates the nth power of a number
real intpow( const real& f, const index& n);
// Does periodic boundary conditions
real pbc( const real& dx, const real& period );

/* REAL VECTOR CLASS */
class rVector {
   friend void pbc(rVector&, const rVector&);
   friend real dotProduct(const rVector&, const rVector&);
   friend std::ostream& operator<<(std::ostream&,const rVector&);
   friend std::istream& operator>>(std::istream&,rVector&);
   friend void chol_elsolve( const rMatrix&, const rVector&, rVector& );
   friend void randNDimGauss( const rMatrix& , const rVector&, rVector& );
   friend void matVecMult( const rMatrix& , const rVector& , rVector& );
   friend void vecMatMult( const rVector& , const rMatrix& , rVector& );
   friend void TmatVecMult( const rMatrix& , const rVector& , rVector& );
   friend void dmatMatMult( const rVector& , const rMatrix& , rMatrix& );
   friend void matDMatMult( const rMatrix&, const rVector& , rMatrix& );
private:
   index sz;
   std::valarray<real> data;
public:
   rVector(index s=0, real dat=0.) : sz(s), data(dat,s) {}
   rVector(const rVector& v) : sz(v.sz), data(v.data) {}
   rVector& operator=(const rVector& v){
      if(&v==this) return *this;
#ifdef DEBUG
      if (v.sz != sz ) ERROR("Vector size mismatch ");
#endif
      data=v.data;
      return *this;
   }
   void resize(const index n){ sz=n; data.resize(sz); } 
   index size() const { return sz; }
   real mag(){ real sum=0.; for(index i=0;i<sz;++i) sum+=data[i]*data[i]; return sqrt(sum); }  
   inline real operator[] (const index& i) const {return data[i];}
   inline real& operator[] (const index& i) {return data[i];}
   rVector& operator=(const real& v) {data=v; return *this;}
   rVector& operator+=(const real& v) {data+=v; return *this;}
   rVector& operator-=(const real& v) {data-=v; return *this;}
   rVector& operator*=(const real& v) {data*=v; return *this;}
   rVector& operator/=(const real& v) {data/=v; return *this;}
   rVector& operator+=(const rVector& v){
#ifdef DEBUG
      if (v.sz != sz ) ERROR("Vector size mismatch ");
#endif
      data+=v.data;
      return *this;
   }
   rVector operator-=(const rVector& v){
#ifdef DEBUG
      if (v.sz != sz ) ERROR("Vector size mismatch ");
#endif
      data-=v.data;
      return *this;
   }
   real max() const {real maxdat=-99999.; for(index i=0;i<sz;++i){ if(data[i]>maxdat) maxdat=data[i]; } return maxdat; }
   real moment(const index& n) const; 
   real central_moment(const index& n) const;
};

// Input and output for rVector
std::ostream& operator<<(std::ostream& ostr,const rVector& vec);
std::istream& operator>>(std::istream& istr,rVector& vec);
// Periodic boundary conditions for a vector
void pbc( rVector& A, const rVector& period );
// Distances between vectors
real calcDistance( const rp::rVector& period, const rp::rVector& v1, const rp::rVector& v2 );
real calcDistance( const rp::rVector& v1, const rp::rVector& v2 );

/* REAL 2D MATRIX CLASS */
class rMatrix {
   friend void matMatMult( const rMatrix& , const rMatrix& , rMatrix& );
   friend void tmatMatMult( const rMatrix& , const rMatrix& , rMatrix& );
   friend void matTMatMult(  const rMatrix& , const rMatrix& , rMatrix& ); 
   friend void dmatMatMult( const rVector& , const rMatrix& , rMatrix& );
   friend void matDMatMult( const rMatrix&, const rVector&, rMatrix& );
#ifndef DRIVER
   friend void invertMat(const rMatrix& , rMatrix& );
   friend void diagMat( const rMatrix&, rVector& , rMatrix&  );
   friend void pseudoInverse( const rMatrix& , rMatrix& );
   friend real logdet( const rMatrix& );
#endif
   friend void transpose( const rMatrix& , rMatrix& );
//   friend real determinant( const rMatrix& );
   friend void cholesky( const rMatrix& , rMatrix& );
   friend void chol_elsolve( const rMatrix&, const rVector&, rVector& );
   friend void randNDimGauss( const rMatrix& , const rVector&, rVector&  );
   friend std::ostream& operator<<(std::ostream& ,const rMatrix& );
   friend std::istream& operator>>(std::istream& ,rMatrix& ); 
   friend void TmatVecMult( const rMatrix& , const rVector& , rVector& );
   friend void matVecMult(const rMatrix& , const rVector& , rVector& );
   friend void vecMatMult( const rVector& , const rMatrix& , rVector& );
private:
   index rw,cl,sz;
   std::valarray<real> data;
public:
   rMatrix(index nr=0, index nc=0, real nv=0.) : rw(nr), cl(nc), sz(nr*nc), data(nv,nr*nc) {}
   rMatrix(const rMatrix& m) : rw(m.rw), cl(m.cl), sz(m.sz), data(m.data) {}
   rMatrix& operator=(const rMatrix& m){ 
      if (&m==this) return *this; 
#ifdef DEBUG
      if (m.rw!=rw || m.cl!=cl) ERROR("Array size mismatch " );
#endif
      data=m.data; 
      return *this; 
   }

   void resize(const index nr=0, const index nc=0) { rw=nr; cl=nc; sz=rw*cl; data.resize(sz); }
   void setToUnity(){ 
#ifdef DEBUG
     if(rw!=cl) ERROR("SQUARE MATRICES CAN'T EQUAL 1"); 
#endif
     data=0.0; for(index i=0;i<rw;++i) data[i+i*cl] = 1.0;
   }
   index rows() const { return rw; }   index cols() const { return cl; }   index size() const { return sz; }
   inline real  operator() (const index& i, const index& j) const { return data[j+i*cl]; }
   inline real& operator() (const index& i, const index& j) { return data[j+i*cl]; }
   real  operator[] (const index& i) const { return data[i]; }
   real& operator[] (const index& i)       { return data[i]; }
   rMatrix& operator=(const real& v) { data=v; return *this;}
   rMatrix& operator+=(const real& v) { data+=v; return *this;}
   rMatrix& operator-=(const real& v) { data-=v; return *this;}
   rMatrix& operator*=(const real& v) { data*=v; return *this;}
   rMatrix& operator/=(const real& v) { data/=v; return *this;}
   rMatrix& operator+=(const rMatrix& m){ 
#ifdef DEBUG
     if (m.rw!=rw || m.cl!=cl) ERROR("Array size mismatch " );
#endif
     data+=m.data;
     return *this; 
   }
   rMatrix& operator-=(const rMatrix& m){ 
#ifdef DEBUG
     if (m.rw!=rw || m.cl!=cl) ERROR("Array size mismatch " );
#endif
     data-=m.data;
     return *this; 
   }
   void setRow(const index& n, const rVector& dat){ 
#ifdef DEBUG
     if(dat.size()!=cl) ERROR("Size mismatch");
#endif
     for(index i=0;i<cl;i++){ data[n*cl+i]=dat[i];} 
   }
   void getRow(const index& n, rVector& dat) const {
#ifdef DEBUG
     if(dat.size()!=cl) ERROR("Size mismatch");
#endif
     for(index i=0;i<cl;i++){ dat[i]=data[n*cl+i];}
   }
};

// Input and output for rMatrix
std::ostream& operator<<(std::ostream& ostr,const rMatrix& mat);
std::istream& operator>>(std::istream& istr,rMatrix& mat);

/* VECTOR VECTOR FUNCTIONS */
// Computes the dot product of two vectors
double dotProduct(const rVector& A, const rVector& B); 
/* MATRIX VECTOR FUNCTIONS */
// Multiplies a vector by a matrix
void matVecMult( const rMatrix& A, const rVector& B, rVector& C);
void vecMatMult( const rVector& A, const rMatrix& B, rVector& C);
// TRANSPOSE(A)*B
void TmatVecMult( const rMatrix& A, const rVector& B, rVector& C);
/* MATRIX FUNCTIONS */
// computes C=A.B matrix matrix multiply
void matMatMult(const rMatrix& A, const rMatrix& B, rMatrix& C);
// Multiply a diagonal matrix by a matrix
void dmatMatMult( const rVector& A, const rMatrix& B, rMatrix& C);
void matDMatMult( const rMatrix& A, const rVector& B, rMatrix& C);
// TRANSPOSE(A)*B
void tmatMatMult(const rMatrix& A, const rMatrix& B, rMatrix& C);
// A*TRANSPOSE(B)
void matTMatMult( const rMatrix& A, const rMatrix& B, rMatrix& C);
// Diagonalizes a symmetric matrix 
void diagMat( const rMatrix& A, rVector& eigenvals, rMatrix& eigenvecs );
// Gets the transpose of a matrix
void transpose(const rMatrix& A, rMatrix& B);
#ifndef DRIVER
// inverts a symmetric matrix
void invertMat(const rMatrix& A, rMatrix& inverse);
// Computes a moore penrose pseudo inverse of a matrix
void pseudoInverse( const rMatrix& A, rMatrix& pseudoinverse );
// Gets the logarithm of the determinant of a matrix
real logdet( const rMatrix& A );
#endif
// Does a cholsky decomposition of a square matrix
void cholesky( const rMatrix& A, rMatrix& B);
void chol_elsolve( const rMatrix& M, const rVector& b, rVector& y );
// Generates a random vector distributed according to a N dimensional gaussian distribution
// with cholesky decomposition M
void randNDimGauss( const rMatrix& M, const rVector& center, rVector& b );
// Generates a Gaussian distributed random number
real randGauss();
};  //ends namespace rp

//!   +++++++++++++++++ STUFF TO CONVERT STRINGS TO INTEGERS AND SO ON ++++++++++++

namespace rp{
    inline double str2float(const std::string& str)
    {
        std::stringstream sstr(str);
        double rval; sstr >> rval;
        return rval;
    }

    inline int str2int(const std::string& str)
    {
        std::stringstream sstr(str);
        int rval; sstr >> rval;
        return rval;
    }

    inline std::string int2str(const long& ival)
    {
        std::stringstream sstr;
        sstr << ival;
        std::string rval; sstr >> rval;
        return rval;
    }

    inline std::string float2str(const double& ival)
    {
        std::stringstream sstr;
        sstr << ival;
        std::string rval; sstr >> rval;
        return rval;
    }
};     // ENDS NAMESPACE RP

//!   +++++++++++++++++ MPI-SAFE I/O UTILITIES +++++++++++++++++++
/************************
  NULL OUTPUT STREAM
*************************/
namespace rp{
struct nullstream: public std::ostream {
    struct nullbuf: std::streambuf {
        int overflow(int c) { return traits_type::not_eof(c); }
    } m_sbuf;
    nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
        };
};

#ifndef _GAT_EXTERN
namespace std{
    extern rp::nullstream cnull;
}; 
#else
namespace std{
    rp::nullstream cnull;
};
#endif

/************************
  MPI SAFE STREAMS
*************************/
namespace rp{    
//MPI-SAFE OUTPUT (execute on every node but output only on node 0)
//also always flushes the data so is best not used when there is intensive i/o activity
class mpiostream: public std::ostream {
private:
    std::ostream& os;
#ifdef MPI
    MPI_Comm mycomm;
#endif
    
public:
#ifdef MPI
    mpiostream(std::ostream& ros, const MPI_Comm& rcomm=MPI_COMM_WORLD) : os(ros), mycomm(rcomm) {}
#else
    mpiostream(std::ostream& ros) : os(ros) {}
#endif
    
    template<class T> std::ostream& operator << (T data) 
    {
#ifdef MPI
        static int myrank=-1;
        if (myrank==-1) MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        if (myrank==0){ os<<data; os.flush(); return os; } else return std::cnull;
#else
        os<<data; os.flush(); return os;
#endif
    }
    
    operator std::ostream&()
    {
#ifdef MPI
        static int myrank=-1;
        if (myrank==-1) MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        if (myrank==0) return os; else return std::cnull;
#else
        return os;
#endif
    }
};
};

#ifndef _GAT_EXTERN
namespace std{
    extern rp::mpiostream pout, perr;
}; 
#endif


#endif //ends #ifndef __rp_utilsh
