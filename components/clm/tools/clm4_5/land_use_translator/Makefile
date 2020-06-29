# first make sure your environment is setup before running 'make' at the command line
# you need a compiler (GNU or INTEL have been tested) and the netcdf library

# Cori
# The INTEL environment is the default
# So all you have to do is: module load cray-netcdf
# Unless you want to switch to the GNU environment: module switch PrgEnv-intel PrgEnv-gnu
# Then run: module load cray-netcdf

# Compy
# No default compiler
# INTEL: module load intel
# GNU: module load gcc
# Then: module load netcdf
# BOTH: make sure that the appropriate LIB_NETCDF and INC_NETCDF are set below for the compiler

# host machine
# check for NERSC_HOST cuz hostname is not unique at NERSC
HN := $(shell hostname)

ifeq ($(NERSC_HOST), cori)
   HN := $(NERSC_HOST)
endif
$(info HN is $(HN))

# optimization flag
OPTIM := -O2

# compiler at NERSC knows where netcdf is because the prg env is set
# but it doesn't know which compiler associated with the current environment
ifeq ($(HN),cori09)

   LIB_NETCDF := $(NETCDF_DIR)/lib

   INC_NETCDF := $(NETCDF_DIR)/include

   # Set environment specified C compiler based on the netcdf compile (gcc or intel)
   ENVCOMP := $(findstring GNU,$(NETCDF_DIR))
   ifeq ($(ENVCOMP),$(NULL))
      ENVCOMP := $(findstring INTEL,$(NETCDF_DIR))
      ifeq ($(ENVCOMP),$(NULL))
         @echo "Use either GNU gcc or INTEL icc compiler"
         exit
      else
         CC_ENV := icc
      endif
   else
      CC_ENV := gcc
   endif
   $(info ENVCOMP is $(ENVCOMP))
endif

# compiler at COMPY does not know where netcdf is
# but it does know which compiler is loaded cuz you have to load it separately
# set this for the default gcc/4.8.5 compiler
# also note the default intel/19.0.3 compiler location
# if you want to use a different compiler you need to change this
ifeq ($(HN),compy-e.pnl.gov)

   LIB_NETCDF := /share/apps/netcdf/4.6.3/gcc/4.8.5/lib
   #LIB_NETCDF := /share/apps/netcdf/4.6.3/intel/19.0.3/lib

   INC_NETCDF := /share/apps/netcdf/4.6.3/gcc/4.8.5/include
   #INC_NETCDF := /share/apps/netcdf/4.6.3/intel/19.0.3/include

   CC_ENV := $(CC)

endif

$(info Compiler (CC_ENV) is $(CC_ENV))

CFLAGS=-lnetcdf -lm -I$(INC_NETCDF) -L$(LIB_NETCDF)

land_use_translator: updateannuallanduse_v2.c
   $(CC_ENV) -std=c11 $(OPTIM) -o land_use_translator updateannuallanduse_v2.c $(CFLAGS)

clean:
	rm -f land_use_translator
