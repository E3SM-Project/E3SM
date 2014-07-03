#-----------------------------------------------------------------------
#
# File:  sgialtix_serial.gnu
#
#  Contains compiler and loader options for the SGI Altix using the 
#  intel compiler and specifies the serial directory for communications 
#  modules.
#
#-----------------------------------------------------------------------
F77 = ifort
F90 = ifort
LD = ifort
CC = cc
Cp = /bin/cp
Cpp = cpp -P
AWK = /usr/bin/gawk
ABI = 
COMMDIR = serial
 
#  Enable MPI library for parallel code, yes/no.

MPI = no

# Adjust these to point to where netcdf is installed

# These have been loaded as a module so no values necessary
#NETCDFINC = -I/netcdf_include_path
#NETCDFLIB = -L/netcdf_library_path
NETCDFINC = -I/usr/projects/climate/maltrud/local/include_8.1
NETCDFLIB = -L/usr/projects/climate/maltrud/local/lib_8.1

#  Enable trapping and traceback of floating point exceptions, yes/no.
#  Note - Requires 'setenv TRAP_FPE "ALL=ABORT,TRACE"' for traceback.

TRAP_FPE = no

#------------------------------------------------------------------
#  precompiler options
#------------------------------------------------------------------

#DCOUPL              = -Dcoupled

Cpp_opts =   \
      $(DCOUPL)

Cpp_opts := $(Cpp_opts) -DPOSIX 
 
#----------------------------------------------------------------------------
#
#                           C Flags
#
#----------------------------------------------------------------------------
 
CFLAGS = $(ABI) 

ifeq ($(OPTIMIZE),yes)
  CFLAGS := $(CFLAGS) -O 
else
  CFLAGS := $(CFLAGS) -g
endif
 
#----------------------------------------------------------------------------
#
#                           FORTRAN Flags
#
#----------------------------------------------------------------------------

FBASE = $(ABI) $(NETCDFINC) -I$(DepDir) -static-libcxa
MODSUF = mod

ifeq ($(TRAP_FPE),yes)
  FBASE := $(FBASE) 
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS = $(FBASE) -O3
else
  FFLAGS = $(FBASE) -g 
endif
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) -v -static-libcxa
 
LIBS = $(NETCDFLIB) -lnetcdf
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------
