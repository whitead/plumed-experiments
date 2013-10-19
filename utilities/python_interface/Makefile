cfiles = $(wildcard *.c)
OBJ= $(cfiles:.c=.o) 


ifeq ($(arch),gnu)
    CC  = gcc -O3 -DAMBER
    LINK = gcc
endif

ifeq ($(arch),dbg)
    CC  = gcc -g -DAMBER
    LINK = gcc
endif

ifeq ($(arch),intel)
    CC  = icc -O3 -DAMBER
    LINK = icc
endif

ifeq ($(arch),'')
  error:
        @echo 'ERROR: Please specify the right architecture (see Makefile for choice)'
endif

# for linux
LDFLAGS= -shared -Wl,-soname,libplumed.so 
# for mac
#LDFLAGS=  -bundle -flat_namespace -undefined suppress   

EXEC = libplumed.so

all:$(EXEC)

gnu:
	@make arch=gnu
dbg:
	@make arch=dbg
intel:
	@make arch=intel

$(EXEC): $(OBJ)
	$(LINK) $(LDFLAGS) $(CFlAGS) $(OBJ) -o $(EXEC)  $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -fPIC -o $@ $*.c

clean:
	@rm -rf *.o  
