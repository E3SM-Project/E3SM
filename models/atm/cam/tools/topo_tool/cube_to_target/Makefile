EXEDIR = .
EXENAME = cube_to_target
RM = rm

.SUFFIXES:
.SUFFIXES: .F90 .o

FC = lf95
DEBUG = FALSE


# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

#------------------------------------------------------------------------
# LF95
#------------------------------------------------------------------------
#
# setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/lib
#
ifeq ($(FC),lf95)
#
# Tramhill
#
  INC_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/include
  LIB_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/lib

  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl -mcmodel=medium 
  FFLAGS   := -c --trace --trap --wide -CcdRR8 -I$(INC_NETCDF)
  ifeq ($(DEBUG),TRUE)
#   FFLAGS += --chk aesu  -Cpp --trace
    FFLAGS += -g --chk a,e,s,u --pca
  else
    FFLAGS += -O
  endif

endif



.F90.o:
	$(FC) $(FFLAGS) $<

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
#reconstruct.o : shr_kind_mod.o
