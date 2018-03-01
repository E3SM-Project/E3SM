#-----------------------------------------------------------------------
#
# File:  linuxabsoft_serial.gnu
#
#  Contains compiler and loader options for the Linux OS using the 
#  Absoft compiler and specifies the serial directory for communications
#  modules.  Use F77 for fortran fixed-form, F90 for free form.
#
#-----------------------------------------------------------------------

F77 = f77
F90 = f95
LD = f95
CC = cc
Cp = /bin/cp
Cpp = cpp
AWK = /usr/bin/gawk
ABI = 
COMMDIR = serial
 
#  Enable MPI library for parallel code, yes/no.

MPI = no

# Adjust these to point to where netcdf is installed

#NETCDFINC = -I/netcdf_include_path
#NETCDFLIB = -L/netcdf_library_path
NETCDFINC = -p/usr/local/include 
NETCDFLIB = -L/usr/local/lib

#  Enable trapping and traceback of floating point exceptions, yes/no.
#  Note - Requires 'setenv TRAP_FPE "ALL=ABORT,TRACE"' for traceback.

TRAP_FPE = no

#------------------------------------------------------------------
#  precompiler options
#------------------------------------------------------------------

#DCOUPL              = -Dcoupled

Cpp_opts =   \
      $(DCOUPL)

Cpp_opts := $(Cpp_opts) -DPOSIX -P -traditional
 
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
 
FBASE = $(ABI) $(NETCDFINC) -p$(DepDir)
MODSUF = mod

ifeq ($(TRAP_FPE),yes)
  FBASE := $(FBASE) 
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS = $(FBASE) -O3
else
  FFLAGS := $(FBASE) -g -Rb -Rc -Rs -trap=ALL
  #FFLAGS := $(FBASE) -g -Rb -Rc -Rs -trap=DIVBYZERO,INVALID
  #FFLAGS := $(FBASE) -g -Rc -Rs -trap=DIVBYZERO,INVALID
  #FFLAGS = $(FBASE) -g 
endif
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------

ifeq ($(OPTIMIZE),yes)
  LDFLAGS = $(ABI) 
else
  LDFLAGS = $(ABI) -debug:full -debug:cv 
endif
 
LIBS = $(NETCDFLIB) -lnetcdf
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
#LDLIBS = $(TARGETLIB) $(LIBRARIES) $(LIBS)
LDLIBS = $(LIBS)
 
