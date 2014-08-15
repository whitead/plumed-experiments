plumed-pinion
=============

A port of PLUMED 1.3.0 for development in the Voth group


CP2K Notes
=========
Using specifically version 13189, I patched the makefile with:

    --- Makefile    (revision 14191)
    +++ Makefile    (working copy)
    @@ -215,16 +215,6 @@
        FCLOGPIPE = &> $*.log
     endif
     
    -#
    -#PLUMED
    -#
    -ifeq ($(PLUMED),yes)
    -       include $(PLUMEDINC) 
    -       DFLAGS += -DPLUMED_CP2K
    -       CPPFLAGS += -DPLUMED_CP2K
    -endif
    -
    -
     # Define the whole bunch of libraries needed
     ALL_LIB =  $(LIB_CP2K_ARCHIVE) $(LIB3_ARCHIVE) $(LIB2_ARCHIVE) $(LIB1_ARCHIVE) $(LIB4_ARCHIVE) $(LIBMA_ARCHIVE) $(LIBNV_ARCHIVE) $(LIB2CUDA_ARCHIVE) $(LIBCUSMM_ARCHIVE)
     ALL_SOURCE_FILES=$(notdir $(OBJECTS:.o=.F) $(LIB1_OBJECTS:.o=.F) $(LIB2_OBJECTS:.o=.F) $(LIB4_OBJECTS:.o=.F) $(LIB3_OBJECTS:.o=.F) cp_common_uses.h cp2k.F cp2k_shell.F)
    @@ -248,7 +238,20 @@
     FCFLAGS += $(CPPSHELL)
     endif
     
    +#
    +#PLUMED
    +#
    +ifeq ($(PLUMED),yes)
    +       include $(PLUMEDINC) 
    +       DFLAGS += -DPLUMED_CP2K
    +       CPPFLAGS += -DPLUMED_CP2K
    +       PLFLAGS:=$(filter-out -traditional,$(CPPFLAGS))
    +       PLFLAGS:=$(filter-out -P,$(PLFLAGS))
    +endif
     
    +
    +
    +
     ### Slave rules ###
     vpath %.F $(SRCDIRS)
     vpath %.h $(SRCDIRS)
    @@ -312,8 +315,9 @@
     # and the rules doing the actual work
     #
     $(OBJ_PLUMED):
    -       $(CC) $(CPPFLAGS) -c $(PLUMEDDIR)/$(@:.o=.cpp)
    +       $(CC) $(PLFLAGS) -c $(PLUMEDDIR)/$(@:.o=.cpp)
     
    +
     _all: $(ALL_LIB) $(PROG) fes $(CP2KSHELL)
     _progr: $(PROG)
     _lib: $(ALL_LIB)


and I used the following arch file:

    #intel/14.0
    #intelmpi/4.1+intel-14.0
    #mkl/11.1
    #
    # export FORT_C_NAME=intel
    # make plumed ARCH=foo PLUMED=yes
    
    XC_LIB = -L${HOME}/xc/lib
    XC_INC = -I${HOME}/xc/include
    
    CC       = icc
    CPP      = cpp
    FC       = mpiifort
    FC_fixed = mpiifort
    LD       = mpiifort
    AR       = ar -r
    DFLAGS   = -D__INTEL -D__FFTSG -D__parallel -D__SCALAPACK -D__BLACS -D__LIBXC2 -D__FFTW3
    CPPFLAGS = -C $(DFLAGS) -P -traditional
    FCFLAGS  = -O2 -mkl -funroll-loops -fpp -free ${XC_INC} ${INT_INC}
    LDFLAGS  = $(FCFLAGS) -L${MKLROOT} ${XC_LIB} ${INT_LIB}
    LIBS     = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -lxc -lstdc++ ${XC_LIB}
    
    OBJECTS_ARCHITECTURE = machine_intel.o
