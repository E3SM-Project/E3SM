#MODEL_FORMULATION = -DNCAR_FORMULATION
MODEL_FORMULATION = -DLANL_FORMULATION

# This flag must be off for nersc hopper:
FILE_OFFSET = -DOFFSET64BIT


dummy:
	( $(MAKE) error )

xlf:
	( $(MAKE) all \
	"FC_PARALLEL = mpxlf90" \
	"CC_PARALLEL = mpcc" \
	"FC_SERIAL = xlf90" \
	"CC_SERIAL = xlc" \
	"FFLAGS_OPT = -O3 -qrealsize=8" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -O0 -g -C -qrealsize=8" \
	"CFLAGS_DEBUG = -O0 -g" \
	"LDFLAGS_DEBUG = -O0 -g" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) $(FILE_OFFSET)" )
 
ftn:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -i4 -r8 -gopt -O2 -Mvect=nosse -Kieee -convert big_endian" \
	"CFLAGS_OPT = -fast" \
	"LDFLAGS_OPT = " \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET)" )

pgi:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = pgf90" \
	"CC_SERIAL = pgcc" \
	"FFLAGS_OPT = -r8 -O3 -byteswapio -Mfree" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -r8 -O0 -g -Mbounds -Mchkptr -byteswapio -Mfree" \
	"CFLAGS_DEBUG = -O0 -g" \
	"LDFLAGS_DEBUG = -O0 -g -Mbounds -Mchkptr" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET)" )

pgi-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -r8 -O3 -byteswapio -Mfree" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET)" )

pgi-llnl:
	( $(MAKE) all \
	"FC_PARALLEL = mpipgf90" \
	"CC_PARALLEL = pgcc" \
	"FC_SERIAL = pgf90" \
	"CC_SERIAL = pgcc" \
	"FFLAGS_OPT = -i4 -r8 -g -O2 -byteswapio" \
	"CFLAGS_OPT = -fast" \
	"LDFLAGS_OPT = " \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET)" )

ifort:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = gcc" \
	"FC_SERIAL = ifort" \
	"CC_SERIAL = gcc" \
	"FFLAGS_OPT = -real-size 64 -O3 -convert big_endian -FR" \
	"CFLAGS_OPT = -O3 -m64" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -real-size 64 -g -convert big_endian -FR -CU -CB -check all" \
	"CFLAGS_DEBUG = -g -m64" \
	"LDFLAGS_DEBUG = -g" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE -m64 $(FILE_OFFSET)" )

gfortran:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = gfortran" \
	"CC_SERIAL = gcc" \
	"FFLAGS_OPT = -O3 -m64 -ffree-line-length-none -fdefault-real-8 -fconvert=big-endian -ffree-form" \
	"CFLAGS_OPT = -O3 -m64" \
	"LDFLAGS_OPT = -O3 -m64" \
	"FFLAGS_DEBUG = -g -m64 -ffree-line-length-none -fdefault-real-8 -fconvert=big-endian -ffree-form -fbounds-check -fbacktrace" \
	"CFLAGS_DEBUG = -g -m64" \
	"LDFLAGS_DEBUG = -g -m64" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE -m64 $(FILE_OFFSET)" )

g95:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = g95" \
	"CC_SERIAL = gcc" \
	"FFLAGS_OPT = -O3 -ffree-line-length-huge -r8 -fendian=big" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET)" )

pathscale-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -r8 -O3 -freeform -extend-source" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET)" )

cray-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -default64 -O3 -f free" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET)" )

intel-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -real-size 64 -O3 -FR" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET)" )

CPPINCLUDES = 
FCINCLUDES = 
LIBS = 
ifneq ($(wildcard $(PIO)/lib), ) # Check for newer PIO version
	CPPINCLUDES = -I../inc -I$(NETCDF)/include -I$(PIO)/include -I$(PNETCDF)/include
	FCINCLUDES = -I../inc -I$(NETCDF)/include -I$(PIO)/include -I$(PNETCDF)/include
	LIBS = -L$(PIO)/lib -L$(PNETCDF)/lib -L$(NETCDF)/lib -lpio -lpnetcdf
else
	CPPINCLUDES = -I../inc -I$(NETCDF)/include -I$(PIO) -I$(PNETCDF)/include
	FCINCLUDES = -I../inc -I$(NETCDF)/include -I$(PIO) -I$(PNETCDF)/include
	LIBS = -L$(PIO) -L$(PNETCDF)/lib -L$(NETCDF)/lib -lpio -lpnetcdf
endif

NCLIB = -lnetcdf
NCLIBF = -lnetcdff
ifneq ($(wildcard $(NETCDF)/lib/libnetcdff.*), ) # CHECK FOR NETCDF4
	LIBS += $(NCLIBF)
endif # CHECK FOR NETCDF4
LIBS += $(NCLIB)

RM = rm -f
CPP = cpp -C -P -traditional
RANLIB = ranlib


ifdef CORE

ifeq "$(DEBUG)" "true"

ifndef FFLAGS_DEBUG
	FFLAGS=$(FFLAGS_OPT)
	CFLAGS=$(CFLAGS_OPT)
	LDFLAGS=$(LDFLAGS_OPT)
	DEBUG_MESSAGE="Debug flags are not defined for this compile group. Defaulting to Optimized flags"
else # FFLAGS_DEBUG IF
	FFLAGS=$(FFLAGS_DEBUG)
	CFLAGS=$(CFLAGS_DEBUG)
	LDFLAGS=$(LDFLAGS_DEBUG)
	override CPPFLAGS += -DMPAS_DEBUG
	DEBUG_MESSAGE="Debugging is on."
endif # FFLAGS_DEBUG IF

else # DEBUG IF
	FFLAGS=$(FFLAGS_OPT)
	CFLAGS=$(CFLAGS_OPT)
	LDFLAGS=$(LDFLAGS_OPT)
	DEBUG_MESSAGE="Debugging is off."
endif # DEBUG IF

ifeq "$(SERIAL)" "true"
	FC=$(FC_SERIAL)
	CC=$(CC_SERIAL)
	SFC=$(FC_SERIAL)
	SCC=$(CC_SERIAL)
	SERIAL_MESSAGE="Serial version is on."
else # SERIAL IF
	FC=$(FC_PARALLEL)
	CC=$(CC_PARALLEL)
	SFC=$(FC_SERIAL)
	SCC=$(CC_SERIAL)
	override CPPFLAGS += -D_MPI
	SERIAL_MESSAGE="Parallel version is on."
endif # SERIAL IF

ifeq "$(USE_PAPI)" "true"
	CPPINCLUDES += -I$(PAPI)/include -D_PAPI
	FCINCLUDES += -I$(PAPI)/include
	LIBS += -L$(PAPI)/lib -lpapi
	PAPI_MESSAGE="Papi libraries are on."
else # USE_PAPI IF
	PAPI_MESSAGE="Papi libraries are off."
endif # USE_PAPI IF

ifeq "$(TAU)" "true"
	LINKER=tau_f90.sh
	CPPINCLUDES += -DMPAS_TAU
	TAU_MESSAGE="TAU Hooks are on."
else
	LINKER=$(FC)
	TAU_MESSAGE="TAU Hooks are off."
endif


####################################################
# Section for adding external libraries and includes
####################################################
ifdef MPAS_EXTERNAL_LIBS
 LIBS += $(MPAS_EXTERNAL_LIBS)
endif
ifdef MPAS_EXTERNAL_INCLUDES
 CPPINCLUDES += $(MPAS_EXTERNAL_INCLUDES)
 FCINCLUDES += $(MPAS_EXTERNAL_INCLUDES)
endif
####################################################


all: mpas_main

mpas_main: 
	cd src; $(MAKE) -j1 FC="$(FC)" \
                 CC="$(CC)" \
                 SFC="$(SFC)" \
                 SCC="$(SCC)" \
                 LINKER="$(LINKER)" \
                 CFLAGS="$(CFLAGS)" \
                 FFLAGS="$(FFLAGS)" \
                 LDFLAGS="$(LDFLAGS)" \
                 RM="$(RM)" \
                 CPP="$(CPP)" \
                 CPPFLAGS="$(CPPFLAGS)" \
                 LIBS="$(LIBS)" \
                 CPPINCLUDES="$(CPPINCLUDES)" \
                 FCINCLUDES="$(FCINCLUDES)" \
                 CORE="$(CORE)"
	if [ ! -e $(CORE)_model ]; then mv src/$(CORE)_model .; fi
	@echo ""
	@echo $(DEBUG_MESSAGE)
	@echo $(SERIAL_MESSAGE)
	@echo $(PAPI_MESSAGE)
	@echo $(TAU_MESSAGE)
clean:
	cd src; $(MAKE) clean RM="$(RM)" CORE="$(CORE)"
	$(RM) $(CORE)_model
error: errmsg

else # CORE IF

all: error
clean: errmsg
error: errmsg
	@echo "************ ERROR ************"
	@echo "No CORE specified. Quitting."
	@echo "************ ERROR ************"
	@echo ""

endif # CORE IF

errmsg:
	@echo ""
	@echo "Usage: $(MAKE) target CORE=[core] [options]"
	@echo ""
	@echo "Example targets:"
	@echo "    ifort"
	@echo "    gfortran"
	@echo "    xlf"
	@echo "    pgi"
	@echo ""
	@echo "Availabe Cores:"
	@cd src; ls -d core_* | grep ".*" | sed "s/core_/    /g"
	@echo ""
	@echo "Available Options:"
#@echo "    SERIAL=true   - builds serial version. Default is parallel version."
	@echo "    DEBUG=true    - builds debug version. Default is optimized version."
	@echo "    USE_PAPI=true - builds version using PAPI for timers. Default is off."
	@echo "    TAU=true      - builds version using TAU hooks for profiling. Default is off."
	@echo ""
	@echo "Ensure that NETCDF, PNETCDF, PIO, and PAPI (if USE_PAPI=true) are environment variables"
	@echo "that point to the absolute paths for the libraries."
	@echo ""

