#ifdef RECONMETAD 
#ifndef __rp_typesh
#define __rp_typesh

// Gromacs specific header files for mpi
#if defined (PLUMED_GROMACS3) || defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45)

   #include "config.h"

   #ifdef PLUMED_GROMACS45
     #ifdef GMX_LIB_MPI
       #define MPI
     #endif
   #else
     #ifdef GMX_MPI
       #define MPI 
     #endif
   #endif
#endif

#define R_PI 3.14159265
#define R_EPS 1.11E-16   // This is machine epsilon for type double

#ifdef MPI
// we DO NOT want C++ bindings!!
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include "mpi.h"
#endif

// This puts an undescore after after the names of lapack routines
#ifndef F77_FUNC
#define F77_FUNC(name,NAME) name ## _
#endif

#define ERROR( MSG ) { std::cerr<<"!ERROR in " << __func__ << ":\n" << MSG <<std::endl; exit(-1); }
#define WARNING( MSG ) { std::cerr << "Warning in " << __func__ << ":\n" << MSG <<std::endl; }

namespace rp {
  typedef unsigned long    index;   //should probably be typedef size_t index to comply with valarrays....
  typedef double           real;
}; //ends namespace rp
#endif //ends #ifndef __rp_typesh
#endif

