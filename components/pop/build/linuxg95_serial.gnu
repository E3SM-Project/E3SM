#-----------------------------------------------------------------------
#
# File:  linuxg95_serial.gnu
#
#  Contains compiler and loader options for the Linux OS using the 
#  open source g95 compiler and specifies the serial directory for 
#  communications modules.
#
#-----------------------------------------------------------------------

F77 = g95
F90 = g95
LD = g95
CC = cc
Cp = /bin/cp
Cpp = /lib/cpp -P
AWK = /usr/bin/gawk
ABI = 
COMMDIR = serial
 
#  Enable MPI library for parallel code, yes/no.

MPI = no

# Adjust these to point to where netcdf is installed

#NETCDFINC = -I/netcdf_include_path
#NETCDFLIB = -L/netcdf_library_path
#NETCDFINC = -I/usr/local/include 
#NETCDFLIB = -L/usr/local/lib
NETCDFINC = -I/home/pwjones/netcdf-3.6.1/include
NETCDFLIB = -L/home/pwjones/netcdf-3.6.1/lib

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
 
FBASE = $(ABI) $(NETCDFINC) -I$(DepDir)
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
 
LDFLAGS = $(ABI) 
 
LIBS = $(NETCDFLIB) -lnetcdf -lf95
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
#LDLIBS = $(TARGETLIB) $(LIBRARIES) $(LIBS)
LDLIBS = $(LIBS)
 
