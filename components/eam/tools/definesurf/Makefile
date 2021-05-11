# Makefile to build definesurf on various platforms
# Note: If netcdf library is not built in the standard location, you must set the environment
# variables INC_NETCDF and LIB_NETCDF

EXEDIR = .
EXENAME = definesurf
RM = rm

.SUFFIXES:
.SUFFIXES: .f90 .o

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

# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

ifeq ($(UNAMEM),CRAY)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
.f90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC = f90
FFLAGS = -64 -c -I$(INC_NETCDF)
LDFLAGS = -64 -L/usr/local/lib64/r4i4 -lnetcdf
.f90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC = f90
FFLAGS = -c -stackvar -f -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
.f90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC = xlf90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
.f90.o:
	$(FC) $(FFLAGS) -qsuffix=f=f90 $<
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
.f90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
ifeq ($(USER_FC),$(null))
FC := pgf90
FFLAGS =  -c -I$(INC_NETCDF) -fast
else
FC := $(USER_FC)
endif
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

ifeq ($(FC),lf95)
FFLAGS = -c --trace --trap -I$(INC_NETCDF) -g
LDFLAGS += -g
endif

.f90.o:
	$(FC) $(FFLAGS) $<
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := ao.o ao_i.o area_ave.o binf2c.o cell_area.o  \
        chkdims.o endrun.o fmain.o handle_error.o inimland.o \
        lininterp.o map_i.o max_ovr.o shr_kind_mod.o sghphis.o sm121.o \
        terrain_filter.o varf2c.o wrap_nf.o interplandm.o map2f.o

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS)  *.mod $(EXEDIR)/$(EXENAME)

ao.o: shr_kind_mod.o
ao_i.o: shr_kind_mod.o
area_ave.o: shr_kind_mod.o
binf2c.o: shr_kind_mod.o
cell_area.o: shr_kind_mod.o
chkdims.o:
endrun.o:
fmain.o: shr_kind_mod.o
handle_error.o:
inimland.o: shr_kind_mod.o
lininterp.o: shr_kind_mod.o
map_i.o: shr_kind_mod.o
max_ovr.o: shr_kind_mod.o
shr_kind_mod.o:
sghphis.o: shr_kind_mod.o
sm121.o: shr_kind_mod.o
terrain_filter.o:
map2f.o:
varf2c.o: shr_kind_mod.o
wrap_nf.o:
interplandm.o: 
