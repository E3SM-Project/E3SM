#
# Makefile for Vertical grid generators
#

ARCH=$(shell uname -s)

.SUFFIXES:
.SUFFIXES: .F $(SUFFIXES)

#####################################################
# Source code vertical grid generation
#  sigma equally spaced
#  hybrid coordinates (18 levels only)
#####################################################

EXE_SIGEQ = sigeqx
SRC_SIGEQ = sigequal.F90
OBJ_SIGEQ = $(SRC_SIGEQ:.F=.o)

EXE_HYB   = hybx
SRC_HYB   = hybgrid.F
OBJ_HYB   = $(SRC_HYB:.F=.o)

###############################################
# Target machines
###############################################
NCOBJ=
ifeq ($(ARCH),SunOS)
  F90  = f90
  FBASE = -I. -I/opt/SUNWhpc/include -e -D$(ARCH) 
  FOPT  = -fast -free -xarch=v8plusa -xchip=ultra2 -dalign -stackvar $(OMPFLAGS)
  FLAGS = $(FOPT) $(FBASE) -C
  LDFLAGS = $(FLAGS)
endif

ifeq ($(ARCH),AIX)
  F90   = xlf90_r
  FBASE = -WF,-D$(ARCH) -qfree=f90 -qdpc=e
  FOPT     = -O3 -qarch=auto -qalias=noaryovrlp 
  FLAGS = $(FOPT) $(FBASE) 
  LDFLAGS = $(FLAGS)

ifneq ('$(NETCDF)','')
  FLAGS += -I/usr/local/lib64/r4i4 -WF,-DNETCDF
  LIBS += -L/usr/local/lib64/r4i4 -lnetcdf
  NCOBJ=netcdf_coord_file.o
endif
endif

ifeq ($(ARCH),Linux)
  F90 = ifort
  FBASE = -D$(ARCH)
  FOPT = -nothreads -O2 -w90 -Vaxlib -free
  FLAGS = $(FOPT) $(FBASE)
  LDFLAGS = $(FLAGS)
endif

ifeq ($(ARCH),Darwin)
  F90     = gfortran
  FBASE   = -D$(ARCH) -ffree-form
  FOPT    = -O3
  FLAGS   = $(FOPT) $(FBASE)
  LDFLAGS = $(FLAGS)
endif


all: $(EXE_SIGEQ) $(EXE_HYB)

.F.o: 
	$(F90) $(FLAGS) -c $<


$(EXE_SIGEQ): $(NCOBJ) $(OBJ_SIGEQ) 
	$(F90) $(LDFLAGS) $(OBJ_SIGEQ) $(NCOBJ) -o $(EXE_SIGEQ) $(LIBS)

$(EXE_HYB): $(NCOBJ) $(OBJ_HYB)
	$(F90) $(LDFLAGS) $(OBJ_HYB) $(NCOBJ) -o $(EXE_HYB) $(LIBS)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.lst
	rm -f $(EXE_HYB) $(EXE_SIGEQ)








