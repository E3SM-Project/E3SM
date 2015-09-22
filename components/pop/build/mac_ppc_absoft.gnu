#----------------------------------------------------------------------
#
#  File:  mac_ppc_absoft.gnu
#
#  Contains compiler and loader options for a Power PC Mac OS using 
#  the Absoft compiler and specifies the serial directory for 
#  communications modules.
#
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#
#  Basic commands - use F77 for fixed form F90 for free form
#
#----------------------------------------------------------------------

F77 = f77
F90 = f95
LD = f95
CC = cc
Cp = /bin/cp
Cpp = /lib/cpp -P
AWK = /usr/bin/awk
ABI = 
COMMDIR = serial
 
#----------------------------------------------------------------------
#
# Set up necessary library and include paths.
#
#----------------------------------------------------------------------

#  netcdf paths

NETCDFINC = -I/usr/local/netcdf/include
NETCDFLIB = -L/usr/local/netcdf/lib

#  Enable MPI library for parallel code, yes/no.

MPI = no

#  Enable trapping and traceback of floating point exceptions, yes/no.

TRAP_FPE = no

#------------------------------------------------------------------
#  Set any precompiler options
#------------------------------------------------------------------

#DCOUPL              = -Dcoupled

Cpp_opts =   \
      $(DCOUPL)

Cpp_opts := $(Cpp_opts) -DPOSIX $(NETCDFINC)
 
#----------------------------------------------------------------------------
#
#                           C Flags
#
#----------------------------------------------------------------------------
 
CFLAGS = 

ifeq ($(OPTIMIZE),yes)
  CFLAGS := $(CFLAGS) 
else
  CFLAGS := $(CFLAGS) -g
endif
 
#----------------------------------------------------------------------------
#
#                           FORTRAN Flags
#
#----------------------------------------------------------------------------
 
FBASE = -v -I/usr/local/absoft/include -p /usr/local/netcdf/mod -p$(DepDir) 
MODSUF = mod

ifeq ($(TRAP_FPE),yes)
  FBASE := $(FBASE) 
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS := $(FBASE) -O2 -cpu:g4 -N11 -round=NEAREST -altiVec
else
  FFLAGS := $(FBASE) -g
endif
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(FFLAGS)
 
LIBS = -L/usr/local/absoft/lib -lnetcdf
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) -lmpi
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) -lfpe
endif
 
#LDLIBS = $(TARGETLIB) $(LIBRARIES) $(LIBS)
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------
