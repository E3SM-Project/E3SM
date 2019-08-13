# Makefile to build interpic on various platforms
# Note: If netcdf library is not built in the standard location, you must set the environment
# variables INC_NETCDF and LIB_NETCDF

EXEDIR = .
EXENAME = interpic
RM = rm

.SUFFIXES:
.SUFFIXES: .f90 .o

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
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC = f90
FFLAGS = -64 -c -trapuv -I$(INC_NETCDF) -g -C -DEBUG:trap_uninitialized=ON
LDFLAGS = -64 -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC = f90
FFLAGS = -c -stackvar -f -I$(INC_NETCDF) -g
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC = xlf90
FFLAGS = -c -I$(INC_NETCDF) -qsuffix=f=f90
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
FC = pgf90
FFLAGS =  -c -I$(INC_NETCDF) -fast
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Darwin
#------------------------------------------------------------------------

ifeq ($(UNAMES),Darwin)
FC = ifort
FFLAGS =  -c -I$(INC_NETCDF) -ftz -g -traceback
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := control.o dimensions.o driver.o err_exit.o fill_positions.o fmain.o \
        handle_special_cases.o interp_driver.o is_special_case.o lininterp.o \
        shr_kind_mod.o varspecs_mod.o wrap_nf.o cpvar.o compare_var.o addglobal.o

.f90.o:
	$(FC) $(FFLAGS) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME)

control.o: control.f90
cpvar.o: cpvar.f90 shr_kind_mod.o
dimensions.o: dimensions.f90 shr_kind_mod.o
driver.o: driver.f90 shr_kind_mod.o varspecs_mod.o control.o dimensions.o fill_positions.o
err_exit.o: err_exit.f90
fill_positions.o: fill_positions.f90 shr_kind_mod.o dimensions.o varspecs_mod.o 
fmain.o: fmain.f90 control.o dimensions.o
handle_special_cases.o: handle_special_cases.f90 shr_kind_mod.o
interp_driver.o: interp_driver.f90 shr_kind_mod.o varspecs_mod.o
is_special_case.o: is_special_case.f90
lininterp.o: lininterp.f90 shr_kind_mod.o
shr_kind_mod.o: shr_kind_mod.f90
varspecs_mod.o: varspecs_mod.f90 shr_kind_mod.o
wrap_nf.o: wrap_nf.f90 shr_kind_mod.o
compare_var.o: compare_var.f90 varspecs_mod.o
addglobal.o: addglobal.f90
