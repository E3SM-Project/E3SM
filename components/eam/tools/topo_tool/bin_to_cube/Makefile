EXEDIR = .
EXENAME = bin_to_cube
RM = rm

.SUFFIXES:
.SUFFIXES: .F90 .o


#
# setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/lib
#

FC = lf95
#DEBUG=TRUE

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

ifeq ($(FC),lf95)
#
# Tramhill
#
  INC_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/include
  LIB_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/lib

  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl -mcmodel=medium 
  FFLAGS   := -c --trace --trap --wide -CcdRR8 -I$(INC_NETCDF)
  ifeq ($(DEBUG),TRUE)
    #TBH:  this works    FFLAGS += -g --chk --pca
    #TBH:  this FAILS    FFLAGS += -g --chk a,e,s,u,x --pca
    FFLAGS += -g --chk a,e,s,u --pca
  else
    FFLAGS += -O
  endif

endif


#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC = xlf90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
.F90.o:
	$(FC) $(FFLAGS) -qsuffix=f=F90 $<
endif


.F90.o:
	$(FC) $(FFLAGS) $<

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := bin_to_cube.o shr_kind_mod.o

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS)  *.mod $(EXEDIR)/$(EXENAME)

bin_to_cube.o: shr_kind_mod.o
