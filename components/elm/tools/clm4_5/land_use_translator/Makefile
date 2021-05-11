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
ifeq ($(HN),cori)

   LIB_NETCDF := $(NETCDF_DIR)/lib

   INC_NETCDF := $(NETCDF_DIR)/include

   # Set environment specified C compiler based on the netcdf compile (gcc or intel)
   ENVCOMP := $(findstring GNU,$(NETCDF_DIR))
   ifeq ($(ENVCOMP),$(NULL))
      ENVCOMP := $(findstring INTEL,$(NETCDF_DIR))
      ifeq ($(ENVCOMP),$(NULL))
         $(info Use either GNU gcc or INTEL icc compiler)
         $(info If you get a make error, ensure that you have used: module load cray-netcdf)
         exit
      else
         CC_ENV := icc
      endif
   else
      CC_ENV := gcc
   endif
   $(info ENVCOMP is $(ENVCOMP))
endif

# compiler at COMPY also knows where netcdf is
# but not which compiler is loaded
ifeq ($(HN),compy-e.pnl.gov)

   LIB_NETCDF := $(NETCDF_ROOT)/lib

   INC_NETCDF := $(NETCDF_ROOT)/include

   # Set environment specified C compiler based on the netcdf compile (gcc or intel)
   ENVCOMP := $(findstring gcc,$(NETCDF_ROOT))
   ifeq ($(ENVCOMP),$(NULL))
      ENVCOMP := $(findstring intel,$(NETCDF_ROOT))
      ifeq ($(ENVCOMP),$(NULL))
         $(info Use either GNU gcc or INTEL icc compiler)
         $(info If you get a make error, ensure that you have used: module load gcc or module load intel, and then module load netcdf)
         exit
      else
         CC_ENV := icc
      endif
   else
      CC_ENV := gcc
   endif
   $(info ENVCOMP is $(ENVCOMP))
endif

$(info Compiler (CC_ENV) is $(CC_ENV))

CFLAGS=-lnetcdf -lm -I$(INC_NETCDF) -L$(LIB_NETCDF)

land_use_translator: updateannuallanduse_v2.c
	$(CC_ENV) -std=c11 $(OPTIM) -o land_use_translator updateannuallanduse_v2.c $(CFLAGS)

clean:
	rm -f land_use_translator
