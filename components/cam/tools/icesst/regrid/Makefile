# Set up special characters
null  :=

# Check for netcdf locations
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

# FFLAGS and LDFLAGS will have at least these values in all situations
FFLAGS  = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

# Determine platform 
UNAMES := $(shell uname -s)

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
  ifeq ($(USER_FC),pgf90)
    FC = pgf90
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g -Mbounds
    else
      FFLAGS += -O
    endif
  else
    FC = lf95
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g --chk esu
    else
      FFLAGS += -O
    endif
  endif
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
  FC = f90
  FFLAGS  += -64 
  LDFLAGS += -64
  ifeq ($(DEBUG),TRUE)
    FFLAGS  += -C -g
    LDFLAGS += -g
  else
    FFLAGS += -O2
  endif
endif

#------------------------------------------------------------------------
# IBM
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
  FC = xlf90
  FFLAGS += -qsuffix=f=f90
  ifeq ($(DEBUG),TRUE)
    FFLAGS  += -C -g
    LDFLAGS += -g
  else
    FFLAGS += -O2
  endif
endif

OBJS = regrid.o binf2c.o wrap_nf.o precision.o map_i.o cell_area.o ao_i.o area_ave.o \
       max_ovr.o ao.o interp_driver.o lininterp.o err_exit.o

regrid: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) regrid *.o *.mod *.stb *.MOD

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) $*.f90

regrid.o:     precision.o
hurrellsst.o: precision.o
binf2c.o:     precision.o
wrap_nf.o:
precision.o:
map_i.o:      precision.o
cell_area.o:  precision.o
ao_i.o:       precision.o
area_ave.o:   precision.o
max_ovr.o:    precision.o
ao.o:         precision.o
interp_driver.o: precision.o
lininterp.o:  precision.o
err_exit.o:
