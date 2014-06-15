# Makefile to build definesurf on various platforms
# Note: If netcdf library is not built in the standard location, you must set the environment
# variables INC_NETCDF and LIB_NETCDF

EXEDIR = .
EXENAME = definehires
RM = rm

.SUFFIXES:
.SUFFIXES: .F90 .o

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

ifeq ($(UNAMEM),CRAY)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
.F90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC = f90
FFLAGS = -c -I$(INC_NETCDF) -64 -mips4 -bytereclen -s -r8 
LDFLAGS = -64 -L/usr/local/lib64/r4i4 -L$(LIB_NETCDF) -lnetcdf
.F90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC = f90
FFLAGS = -c -stackvar -f -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
.F90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC = xlf90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
.F90.o:
	$(FC) $(FFLAGS) -qsuffix=f=F90 $<
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
.F90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
ifeq ($(USER_FC),$(null))
FC := pgf90
FFLAGS =  -c -I$(INC_NETCDF) -fast -r8 -byteswapio
else
FC := $(USER_FC)
endif
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

ifeq ($(FC),lf95)
  FFLAGS   := -c --trace --trap --wide -CcdRR8 -I$(INC_NETCDF)
  ifeq ($(DEBUG),TRUE)
    #TBH:  this works    FFLAGS += -g --chk --pca
    #TBH:  this FAILS    FFLAGS += -g --chk a,e,s,u,x --pca
    FFLAGS += -g --chk a,e,s,u --pca
  else
    FFLAGS += -O
  endif
endif

.F90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := gtopo30_to_10min.o shr_kind_mod.o

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS)  *.mod $(EXEDIR)/$(EXENAME)

gtopo30_to_10min.o: shr_kind_mod.o
