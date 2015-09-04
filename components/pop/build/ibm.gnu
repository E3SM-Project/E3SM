#-----------------------------------------------------------------------
#
# File:  ibm.gnu
#
#  Contains compiler and loader options for the ibm using the 
#  xlf90 compiler and specifies the serial directory for 
#  communications modules.
#
#-----------------------------------------------------------------------

F77 = xlf90_r
F90 = xlf90_r
LD = xlf90_r
CC = cc_r
Cp = /bin/cp
Cpp = /usr/bin/cpp -P
AWK = /usr/bin/awk
ABI = 
COMMDIR = serial
 
#  Enable MPI library for parallel code, yes/no.

MPI = no

# Adjust these to point to where netcdf is installed

NETCDFINC = -I/usr/local/include
NETCDFLIB = /usr/local/lib

#  Enable trapping and traceback of floating point exceptions, yes/no.
#  Note - Requires 'setenv TRAP_FPE "ALL=ABORT,TRACE"' for traceback.

TRAP_FPE = no

#------------------------------------------------------------------
#  precompiler options
#------------------------------------------------------------------

DCOUPL              = -Dcoupled

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
 
FBASE = $(ABI) $(NETCDFINC) -I$(DepDir)
MODSUF = mod

ifeq ($(TRAP_FPE),yes)
  FBASE := $(FBASE) 
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS = $(FBASE) -O2
else
  FFLAGS = $(FBASE) -g 
endif
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = -q64 -bdatapsize:64K -bstackpsize:64K -btextpsize:64K
 
LIBS = -L/usr/local/lib -l/usr/local/lib 
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
#LDLIBS = $(TARGETLIB) $(LIBRARIES) $(LIBS)
LDLIBS = -L/usr/local/lib
 
