# Location of RTE+RRTMGP libraries, module files.
RRTMGP_DIR = ../../build
# Sets macros FC, FCFLAGS consistent with RTE+RRTMGP
include $(RRTMGP_DIR)/Makefile.conf
#
# RRTMGP library, module files
#
LDFLAGS   += -L$(RRTMGP_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_DIR)

#
# netcdf library, module files
# C and Fortran interfaces respectively
#
NCHOME = /opt/local
NFHOME = $(HOME)/Applications/$(FC)
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LIBS      += -lnetcdff -lnetcdf

#
# Setting macro TIMING=yes uses routines from the General Purpose Timing Library
#  https://jmrosinski.github.io/GPTL/
#
TIMING=no
# Compiler specific - e.g. turn off for pgfortran, set to -cpp for gfortran
FCFLAGS += -fpp
ifeq ($(TIMING),yes)
	#
	# Timing library
	#
	FCINCLUDE += -I$(TIME_DIR)/include
	# Compiler specific
	FCFLAGS += -DUSE_TIMING
	# Installation specific
	TIME_DIR = $(HOME)/Codes/GPTL/gptl-v5.5.3/macos
	LDFLAGS   += -L$(TIME_DIR)/lib
	LIBS      += -lgptl
endif

VPATH = ../

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

#
# Ancillary codes
#
ADDITIONS = mo_simple_netcdf.o mo_rfmip_io.o mo_load_coefficients.o

all: rrtmgp_rfmip_lw rrtmgp_rfmip_sw

rrtmgp_rfmip_lw:   rrtmgp_rfmip_lw.o   $(ADDITIONS)

rrtmgp_rfmip_lw.o: rrtmgp_rfmip_lw.F90 $(ADDITIONS)

rrtmgp_rfmip_sw:   rrtmgp_rfmip_sw.o   $(ADDITIONS)

rrtmgp_rfmip_sw.o: rrtmgp_rfmip_sw.F90 $(ADDITIONS)

mo_rfmip_io.o:          mo_rfmip_io.F90          mo_simple_netcdf.o

mo_load_coefficients.o: mo_load_coefficients.F90 mo_simple_netcdf.o

clean:
	-rm rrtmgp_rfmip_lw *.o *.mod *.optrpt ../*.optrpt
