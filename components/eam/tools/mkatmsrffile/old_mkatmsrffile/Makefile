#-------------------------------------
# Makefile for mkatmsrffile tool 
#-------------------------------------

#-------------------------------------
#Instructions:
# Load a compiler and the associated 
# netcdf library. Update variables "FC"
# and "NETCDF_PATH" to add the compiler
# and path to netcdf home directory
#-------------------------------------

ifeq ($(FC),$(null))
   FC := ifort
endif

NETCDF_PATH := $(NETCDF_HOME)

EXE := mkatmsrffile

#---------------------------
# USER INPUT ENDS HERE
#---------------------------

#include string
INC = -I$(NETCDF_ROOT)/include
LIB = -L$(NETCDF_ROOT)/lib -lnetcdf -lnetcdff

FFLAGS += -O2

$(EXE): mkatmsrffile.F90 
	$(FC) $(FFLAGS) $(INC)  $< -o $@  $(LIB)

clean:
	/bin/rm $(EXE)

