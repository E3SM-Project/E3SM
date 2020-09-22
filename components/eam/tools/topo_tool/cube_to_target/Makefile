EXEDIR = .
EXENAME = cube_to_target
RM = rm

.SUFFIXES:
.SUFFIXES: .F90 .o

# Set the compiler; note that on some platforms, FC may be defined as an
# environment variable, and we should not need to do anything here. If this
# is the case, skip setting FC, but if FC is null then set to a default.
ifeq ($(FC),$(null))
    FC = gfortran
endif
DEBUG = FALSE

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
    LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
    INC_NETCDF := /usr/local/include
endif

# Set default compile and link flags
FFLAGS  = -L$(LIB_NETCDF) -I$(INC_NETCDF)
#LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl -mcmodel=medium 
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff -lcurl

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

#------------------------------------------------------------------------
# Special flags for LF95
#------------------------------------------------------------------------
#
ifeq ($(FC),lf95)
  #
  # Tramhill
  #
  INC_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/include
  LIB_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/lib

  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl
-mcmodel=medium 
  FFLAGS   := -c --trace --trap --wide -CcdRR8 -I$(INC_NETCDF)
  ifeq ($(DEBUG),TRUE)
    #FFLAGS += --chk aesu  -Cpp --trace
    FFLAGS += -g --chk a,e,s,u --pca
  else
    FFLAGS += -O
  endif

endif

.F90.o:
	$(FC) -c $(FFLAGS) $<

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := reconstruct.o remap.o cube_to_target.o shr_kind_mod.o

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS)  *.mod $(EXEDIR)/$(EXENAME)

cube_to_target.o: shr_kind_mod.o remap.o reconstruct.o
remap.o: 
reconstruct.o: remap.o
