ffiles = $(wildcard *.f90)
cfiles = $(wildcard *.c)
cxxfiles = $(wildcard *.cpp)
OBJ= $(ffiles:.f90=.o) $(cfiles:.c=.o) $(cxxfiles:.cpp=.o)

arch=g95


ifeq ($(arch),gfortran)
    F90 = gfortran -O3 -fno-second-underscore -g
    CC  = gcc -O3 -DDRIVER -g
    CXX = g++ -O3 -DDRIVER -g -DDEBUG
    LINK = gfortran 
    LIBS = 
endif

ifeq ($(arch),g95)
    F90 = g95 -g -fno-second-underscore     
    CC  = gcc -g -DDRIVER 
    CXX = g++ -g -DDRIVER 
    LINK = g95
    LIBS = 
endif

ifeq ($(arch),gfortran_64)
    F90 = gfortran -m64 -O3 -fno-second-underscore     
    CC  = gcc -O3 -DDRIVER -m64
    CXX = g++ -O3 -DDRIVER 
    LINK = gfortran -m64
    LIBS = -lstdc++ 
endif

ifeq ($(arch),intel)
    F90 = ifort -O3 
    CC  = icc -O3 -DDRIVER 
    CXX = icpc -O3 -DDRIVER 
    LINK = ifort 
    LIBS = 
endif

ifeq ($(arch),'')
  error:
        @echo 'ERROR: Please specify the right architecture (see Makefile for choice)'
endif


EXEC = driver 

all:$(EXEC)

g95:
	@make arch=g95

gfortran:
	@make arch=gfortran

gfortran_64:
	@make arch=gfortran_64

intel:
	@make arch=intel

$(EXEC): $(OBJ)
	$(LINK) $(CFlAGS) $(OBJ) -o $(EXEC)  $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c  -o $@ $*.c

%.o: %.f90
	$(F90) $(PFLAGS) -c  -o $@ $^

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $*.cpp

clean:
	@rm -rf *.o  
