#-----------------------------------------------------------------------
# This Makefile is for building cprnc on AIX, Compaq, Linux (with pgf90,
# lf95, ifort, or gfortran compilers), IRIX or SUN platforms.
#
# These macros can be changed by setting environment variables:
#
# Set the path to netcdf:
#
# gmake NETCDF=pathToNetcdf
#
# This sets LIB_NETCDF=$NETCDF/lib and INC_NETCDF=$NETCDF/include
#
# LIB_NETCDF --- Library directory location of netcdf.
# INC_NETCDF --- Include directory location of netcdf.  This directory needs to contain
#                the NetCDF .mod files, i.e., netcdf.mod and typesizes.mod.
#
# If the include and library files don't have a common root directory then set them
# independently in the commandline:
#
# gmake LIB_NETCDF=pathToLIBFiles INC_NETCDF=pathToINCFiles
#
# You also can set the environment variables:
#
# USER_FC ------ User defined Fortran compiler (for Linux can be pgf90, lf95, ifort, gfortran)
# EXEDIR ------- Directory to build executable in. (Defaults to .)
# VPATH -------- GNU make path. (Defaults to current directory)
#
#------------------------------------------------------------------------
include Macros.make
# Set up special characters
null  :=

EXENAME = cprnc
RM = rm

NETCDF_PATH ?= $(NETCDF_FORTRAN_PATH)

# Default for the netcdf library and include directories
LIB_NETCDF := $(NETCDF_PATH)/lib
INC_NETCDF := $(NETCDF_PATH)/include
LDFLAGS = -L$(LIB_NETCDF) -lnetcdff -Wl,-rpath $(LIB_NETCDF)

# Determine platform
UNAMES := $(shell uname -s)
SNAME := $(shell uname -n | cut -c1-2)

GENF90 = ../../src/externals/genf90/genf90.pl

FC := $(SFC)
FFLAGS += -I$(INC_NETCDF) -I.
#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

# If path to source code not given
ifeq ($(VPATH),$(null))
  VPATH:= .
endif

OBJS := compare_vars_mod.o cprnc.o filestruct.o prec.o utils.o

# If executable directory not given
ifeq ($(EXEDIR),$(null))
  EXEDIR  := .
endif

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o .in

.F90.o:
	$(FC) -c $(FFLAGS) $<

.f90.o:
	$(FC) -c $(FFLAGS) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS) $(SLIBS)

compare_vars_mod.F90 : compare_vars_mod.F90.in
	perl $(GENF90) $< > $@

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME)

# remove generated file during clean
realclean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME) compare_vars_mod.F90 core

include $(CURDIR)/Depends

# 'make check' will run the standard tests but without baseline
# comparisons. For complete testing, you should generally also do
# baseline comparisons. See the notes in test_inputs/README for details.
check: $(EXEDIR)/$(EXENAME)
	$(RM) -fr tmpdir
	./run_tests -outdir tmpdir
