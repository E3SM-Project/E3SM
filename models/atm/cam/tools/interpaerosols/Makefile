# Makefile to build interpaerosols on various platforms
# Note: If netcdf library is not built in the standard location, you must set the environment
# variables INC_NETCDF and LIB_NETCDF
# This program now looks for the file netcdf.mod instead of netcdf.inc
# on some NCAR systems this means that INC_NETCDF and LIB_NETCDF will be set to 
# the same directory
# Set up special characters
null  :=

EXEDIR = .
EXENAME = interpaerosols
RM = rm

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

# Determine platform 
UNAMES := $(shell uname -s)

# Architecture-specific flags and rules
#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC      = f90
FFLAGS  = -64 -c -I$(INC_NETCDF) -O2 -D$(UNAMES) -DHIDE_MPI
LDFLAGS = -64 -L$(LIB_NETCDF) -lnetcdf -lscs
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC      = xlf90
FFLAGS  = -c -I$(INC_NETCDF) -q64 -qsuffix=f=f90:cpp=F90 -O2 -qmaxmem=-1 -g -qfullpath -WF,-DHIDE_MPI,-DAIX
LDFLAGS = -L$(LIB_NETCDF) -q64 -lnetcdf -llapack -lblas -g
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC      = f90
FFLAGS  = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lcxml
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
  ifeq ($(USER_FC),$(null))
    FC := pgf90
  else
    FC := $(USER_FC)
  endif

  FFLAGS = -c -I$(INC_NETCDF) -DHIDE_MPI -D$(UNAMES)
  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

  ifeq ($(FC),pgf90)
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g -Ktrap=fp -Mbounds
    else
      FFLAGS += -fast
    endif
    LDFLAGS += -llapack -lblas
  endif

  ifeq ($(FC),lf95)
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g --chk e,s,u
    else
      FFLAGS += -O
    endif
    LDFLAGS += -llapackmt -lblasmt
  endif
endif

#------------------------------------------------------------------------
# Cray X1
#------------------------------------------------------------------------

ifeq ($(UNAMES),UNICOS/mp)
FC      = ftn
FFLAGS  = -c -I$(INC_NETCDF) -O2 -DUNICOSMP
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := addglobal.o driver.o fmain.o globals.o preserve_mean.o error_messages.o\
        interpolate_data.o shr_kind_mod.o abortutils.o shr_sys_mod.o shr_mpi_mod.o 

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

.f90.o:
	$(FC) $(FFLAGS) $<

.F90.o:
	$(FC) $(FFLAGS) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME)

addglobal.o: error_messages.o
interpolate_data.o: abortutils.o shr_kind_mod.o
driver.o: shr_kind_mod.o globals.o preserve_mean.o interpolate_data.o
fmain.o: globals.o
preserve_mean.o: shr_kind_mod.o globals.o
abortutils.o: shr_sys_mod.o
shr_sys_mod.o: shr_mpi_mod.o
error_messages.o: abortutils.o
shr_mpi_mod.o: shr_kind_mod.o

VPATH = . ../../src/control ../../../../csm_share/shr ../../src/utils
