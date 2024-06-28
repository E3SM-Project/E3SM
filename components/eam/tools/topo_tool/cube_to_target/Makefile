EXEDIR = .
EXENAME = cube_to_target
RM = rm

.SUFFIXES:
.SUFFIXES: .F90 .o

# Set the compiler
FC := gfortran

# Set NetCDF library and include directories
LIB_NETCDF := $(shell nf-config --prefix)/lib
INC_NETCDF := $(shell nf-config --prefix)/include

# Set default compile and link flags
FFLAGS  += -L$(LIB_NETCDF) -I$(INC_NETCDF)
LDFLAGS += -L$(LIB_NETCDF) -lnetcdf -lnetcdff

ifeq ($(FC),gfortran)
    FFLAGS += -fallow-argument-mismatch
endif

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

#------------------------------------------------------------------------
#------------------------------------------------------------------------

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
