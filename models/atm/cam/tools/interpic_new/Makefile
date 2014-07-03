# Makefile to build interpic on various platforms
# Note: If netcdf library is not built in the standard location, you must set the
# variables INC_NETCDF and LIB_NETCDF

EXEDIR = .
EXENAME = interpic
RM = rm

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 

# Check whether NetCDF library and include directories specified in environment
# or on make commandline.
ifeq ($(strip $(LIB_NETCDF)),)
  LIB_NETCDF := /usr/local/lib
endif
ifeq ($(strip $(INC_NETCDF)),)
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
FFLAGS = -c -I$(INC_NETCDF) -qsuffix=cpp=F90 -WF,-DAIX -g -qfullpath
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

# g95
#FC = g95 
#FFLAGS =  -c -I$(INC_NETCDF) -g -ftrace=full

# pgf90
#FC = pgf90
#FFLAGS =  -c -I$(INC_NETCDF) -g -Ktrap=fp -Mrecursive -Mbounds

# lf95
#FC = lf95
#FFLAGS =  -c -I$(INC_NETCDF) -g --chk a,e,s,u --pca --trace --trap

# ifort
FC = ifort
FFLAGS =  -c -I$(INC_NETCDF) -g -check all -fpe0 -traceback -ftz -convert big_endian -fp-model precise

LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := addglobal.o compare_var.o control.o cpvar.o dimensions.o driver.o \
        err_exit.o fill_positions.o fmain.o \
        handle_special_cases.o interp.o interpolate_data.o is_special_case.o \
        shr_kind_mod.o wrap_nf.o

.F90.o:
	$(FC) $(FFLAGS) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME)

addglobal.o: addglobal.F90
compare_var.o: compare_var.F90 fill_positions.o
control.o: control.F90
cpvar.o: cpvar.F90 shr_kind_mod.o
dimensions.o: dimensions.F90 shr_kind_mod.o control.o
driver.o: driver.F90 shr_kind_mod.o control.o dimensions.o fill_positions.o interp.o
err_exit.o: err_exit.F90
fill_positions.o: fill_positions.F90 shr_kind_mod.o control.o
fmain.o: fmain.F90 control.o dimensions.o
handle_special_cases.o: handle_special_cases.F90 shr_kind_mod.o
interp.o: interp.F90 shr_kind_mod.o fill_positions.o interpolate_data.o control.o
interpolate_data.o: interpolate_data.F90
is_special_case.o: is_special_case.F90
shr_kind_mod.o: shr_kind_mod.F90
wrap_nf.o : wrap_nf.F90
