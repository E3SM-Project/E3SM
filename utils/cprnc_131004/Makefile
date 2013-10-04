#-----------------------------------------------------------------------
# This Makefile is for building cprnc on AIX, Compaq, Linux (with pgf90 or 
# lf95 compiler), IRIX or SUN platforms.
#
# These macros can be changed by setting environment variables:
#
# Set the path to netcdf:
# 
# gmake NETCDF=pathToNetcdf
#
# This sets LIB_NETCDF=$NETCDF/lib and INC_NETCDF=$NETCDF/include
#
# LIB_NETCDF --- Library directory location of netcdf. (defaults to /usr/local/lib)
# INC_NETCDF --- Include directory location of netcdf. (defaults to /usr/local/include)
#
# You also can set the environment variables:
#
# USER_FC ------ User defined Fortran compiler (for Linux can be pgf90, lf95, or ifort)
# EXEDIR ------- Directory to build executable in. (Defaults to .)
# VPATH -------- GNU make path. (Defaults to current directory)
#
#------------------------------------------------------------------------

# Set up special characters
null  :=

EXENAME = cprnc
RM = rm

ifeq ($(NETCDF),$(null))
NETCDF := /contrib/netcdf/4.1.3_seq
endif

# Check for the netcdf library and include directories 
LIB_NETCDF := $(NETCDF)/lib
INC_NETCDF := $(NETCDF)/include

# Determine platform 
UNAMES := $(shell uname -s)
SNAME := $(shell uname -n | cut -c1-2)

GENF90 = ./genf90/genf90.pl

# Architecture-specific flags and rules
#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC      = f90
FFLAGS  = -64 -c -r8 -i4 -I$(INC_NETCDF) -O2
LDFLAGS = -64 -L$(LIB_NETCDF) -lnetcdf 
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC      = xlf90
FFLAGS  = -c -g -I$(INC_NETCDF) -q64 -qsuffix=f=f90:cpp=F90 -O2 -qmaxmem=-1
LDFLAGS = -L$(LIB_NETCDF) -q64 -lnetcdf 
endif

#------------------------------------------------------------------------
# Darwin
#------------------------------------------------------------------------

ifeq ($(UNAMES),Darwin)
FC      = xlf90
FFLAGS  = -c -I$(INC_NETCDF) -qsuffix=f=f90 -O2 -qmaxmem=-1
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Cray X1E with UNICOS/mp
#------------------------------------------------------------------------

ifeq ($(UNAMES),UNICOS/mp)
FC      = ftn
FFLAGS  = -c -DUNICOSMP -DCPP_VECTOR -f free -F -d y -I$(INC_NETCDF) -I. -p. -Ovector3,stream3,scalar3
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC      = f90
FFLAGS  = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#-----------------------------------------------------------------------
# SUN
#-----------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC      = f90
FFLAGS  = -c -I$(INC_NETCDF) -stackvar
LDFLAGS = -L$(LIB_NETCDF) -L$(HOME)/lib -lnetcdf
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

  ifeq ($(FC),pgf90)
    ifeq ($(DEBUG),TRUE)
      FFLAGS = -c -I$(INC_NETCDF) -g -Ktrap=fp -Mbounds
    else
      FFLAGS = -c -I$(INC_NETCDF) -fast
    endif
  endif

  ifeq ($(FC),lf95)
    ifeq ($(DEBUG),TRUE)
# currently lf95 debug will fail due to netcdf padding characters
      FFLAGS =  -c -I$(INC_NETCDF) -g --chk a,e,s,u
    else
      FFLAGS =  -c -I$(INC_NETCDF) -O
    endif

  endif

  ifeq ($(FC),ifort)
    FFLAGS =  -c -I$(INC_NETCDF) -132 -ftz -g -m64
    ifeq ($(DEBUG),TRUE)
      FFLAGS +=  -CB
    else
      FFLAGS +=  -O2
    endif
  endif
  ifeq ($(FC),gfortran)
    FFLAGS = -c -I$(INC_NETCDF) -O -ffree-form -ffree-line-length-none
  endif
  LDFLAGS = -L$(LIB_NETCDF) -lnetcdff -lnetcdf


endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

# If path to source code not given
ifeq ($(VPATH),$(null))
  VPATH:= .
endif
dirs := $(subst :,$(space),$(VPATH))

# Get list of files and determine objects and dependency files
FIND_FILES = $(wildcard $(dir)/*.[Ff]90)
FILES      = $(foreach dir, $(dirs),$(FIND_FILES))
SOURCES   := $(sort $(notdir $(FILES)))
OBJS      := $(addsuffix .o, $(basename $(SOURCES)))

ifeq ($(SNAME),be)
	# Need to add this object since it is code generated 
	# won't be in the inital obj list.
	LDFLAGS += -lnetcdff
	OBJS += compare_vars_mod.o
endif
ifeq ($(SNAME),ly)
	OBJS += compare_vars_mod.o
endif

# If executable directory not given
ifeq ($(EXEDIR),$(null))
  EXEDIR  := .
endif

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o .in

.F90.o:
	$(FC) $(FFLAGS) $<

.f90.o:
	$(FC) $(FFLAGS) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

compare_vars_mod.F90 : compare_vars_mod.F90.in
	perl $(GENF90) $< > $@

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME) 
	 
# remove generated file during clean
realclean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME) compare_vars_mod.F90 core

include $(CURDIR)/Depends

