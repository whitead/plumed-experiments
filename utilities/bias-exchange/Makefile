#
# Makefile for the compilation of exchange-tool.x
# change the fortran compiler if necessary
#

FC=gfortran

all: exchange-tool.x

exchange-tool.x: exchange-tool.f90
	$(FC) $(FFLAGS) exchange-tool.f90 -o exchange-tool.x

clean: 
	rm -rf exchange-tool.x
