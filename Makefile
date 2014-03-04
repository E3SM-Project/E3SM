#MODEL_FORMULATION = -DNCAR_FORMULATION
MODEL_FORMULATION = -DLANL_FORMULATION


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
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI" )
 
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
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

pgi:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = pgf90" \
	"CC_SERIAL = pgcc" \
	"FFLAGS_OPT = -r8 -O3 -byteswapio -Mfree" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -r8 -O0 -g -Mbounds -Mchkptr -byteswapio -Mfree -Ktrap=divz,fp,inv,ovf -traceback" \
	"CFLAGS_DEBUG = -O0 -g -traceback" \
	"LDFLAGS_DEBUG = -O0 -g -Mbounds -Mchkptr -Ktrap=divz,fp,inv,ovf -traceback" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

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
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

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
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

ifort:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = ifort" \
	"CC_SERIAL = icc" \
	"FFLAGS_OPT = -real-size 64 -O3 -convert big_endian -FR" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -real-size 64 -g -convert big_endian -FR -CU -CB -check all -fpe0 -traceback" \
	"CFLAGS_DEBUG = -g -fpe0 -traceback" \
	"LDFLAGS_DEBUG = -g -fpe0 -traceback" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

ifort-gcc:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = ifort" \
	"CC_SERIAL = gcc" \
	"FFLAGS_OPT = -real-size 64 -O3 -convert big_endian -FR" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -real-size 64 -g -convert big_endian -FR -CU -CB -check all -fpe0 -traceback" \
	"CFLAGS_DEBUG = -g -traceback" \
	"LDFLAGS_DEBUG = -g -fpe0 -traceback" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

gfortran:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = gfortran" \
	"CC_SERIAL = gcc" \
	"FFLAGS_OPT = -O3 -m64 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fconvert=big-endian -ffree-form" \
	"CFLAGS_OPT = -O3 -m64" \
	"LDFLAGS_OPT = -O3 -m64" \
	"FFLAGS_DEBUG = -g -m64 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fconvert=big-endian -ffree-form -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow" \
	"CFLAGS_DEBUG = -g -m64" \
	"LDFLAGS_DEBUG = -g -m64" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

gfortran-openmpi:
	( $(MAKE) all \
	"FC_PARALLEL = openmpif90" \
	"CC_PARALLEL = openmpicc" \
	"FC_SERIAL = gfortran" \
	"CC_SERIAL = gcc" \
	"FFLAGS_OPT = -O3 -m64 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fconvert=big-endian -ffree-form" \
	"CFLAGS_OPT = -O3 -m64" \
	"LDFLAGS_OPT = -O3 -m64" \
	"FFLAGS_DEBUG = -g -m64 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fconvert=big-endian -ffree-form -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow" \
	"CFLAGS_DEBUG = -g -m64" \
	"LDFLAGS_DEBUG = -g -m64" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

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
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

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
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

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
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

intel-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -real-size 64 -O3 -convert big_endian -FR" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

bluegene:
	( $(MAKE) all \
	"FC_PARALLEL = mpixlf95_r" \
	"CC_PARALLEL = mpixlc_r" \
	"FC_SERIAL = bgxlf95_r" \
	"CC_SERIAL = bgxlc_r" \
	"FFLAGS_OPT = -O2 -g -qrealsize=8" \
	"CFLAGS_OPT = -O2 -g" \
	"LDFLAGS_OPT = -O2 -g" \
	"FFLAGS_DEBUG = -O0 -g -C -qinitalloc -qinitauto -qrealsize=8" \
	"CFLAGS_DEBUG = -O0 -g" \
	"LDFLAGS_DEBUG = -O0 -g" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI" )

CPPINCLUDES = 
FCINCLUDES = 
LIBS = 
ifneq ($(wildcard $(PIO)/lib), ) # Check for newer PIO version
	CPPINCLUDES = -I$(NETCDF)/include -I$(PIO)/include -I$(PNETCDF)/include
	FCINCLUDES = -I$(NETCDF)/include -I$(PIO)/include -I$(PNETCDF)/include
	LIBS = -L$(PIO)/lib -L$(PNETCDF)/lib -L$(NETCDF)/lib -lpio -lpnetcdf
else
	CPPINCLUDES = -I$(NETCDF)/include -I$(PIO) -I$(PNETCDF)/include
	FCINCLUDES = -I$(NETCDF)/include -I$(PIO) -I$(PNETCDF)/include
	LIBS = -L$(PIO) -L$(PNETCDF)/lib -L$(NETCDF)/lib -lpio -lpnetcdf
endif

NCLIB = -lnetcdf
NCLIBF = -lnetcdff
ifneq ($(wildcard $(NETCDF)/lib/libnetcdff.*), ) # CHECK FOR NETCDF4
	LIBS += $(NCLIBF)
endif # CHECK FOR NETCDF4
LIBS += $(NCLIB)

RM = rm -f
CPP = cpp -P -traditional
RANLIB = ranlib

ifdef CORE

ifneq ($(wildcard src/core_$(CORE)), ) # CHECK FOR EXISTENCE OF CORE DIRECTORY

ifneq ($(wildcard src/core_$(CORE)/build_options.mk), ) # Check for build_options.mk
include src/core_$(CORE)/build_options.mk
else # ELSE Use Default Options
EXE_NAME=$(CORE)_model
NAMELIST_SUFFIX=$(CORE)
endif

override CPPFLAGS += -DMPAS_NAMELIST_SUFFIX=$(NAMELIST_SUFFIX)
override CPPFLAGS += -DMPAS_EXE_NAME=$(EXE_NAME)

else # ELSE CORE DIRECTORY CHECK

report_builds: all

endif # END CORE DIRECTORY CHECK

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

FC=$(FC_PARALLEL)
CC=$(CC_PARALLEL)
SFC=$(FC_SERIAL)
SCC=$(CC_SERIAL)
PARALLEL_MESSAGE="Parallel version is on."

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

ifeq "$(GEN_F90)" "true"
	GEN_F90_MESSAGE="MPAS generated and was built with intermediate .f90 files."
else
	override GEN_F90=false
	GEN_F90_MESSAGE="MPAS was built with .F files."
endif


ifneq ($(wildcard .mpas_core_*), ) # CHECK FOR BUILT CORE

ifneq ($(wildcard .mpas_core_$(CORE)), ) # CHECK FOR SAME CORE AS ATTEMPTED BUILD.
	override AUTOCLEAN=false
	CONTINUE=true
else
	LAST_CORE=`cat .mpas_core_*`

ifeq "$(AUTOCLEAN)" "true" # CHECK FOR CLEAN PRIOR TO BUILD OF A NEW CORE.
	CONTINUE=true
	AUTOCLEAN_MESSAGE="Infrastructure was cleaned prior to building ."
else
	CONTINUE=false
endif # END OF AUTOCLEAN CHECK

endif # END OF CORE=LAST_CORE CHECK

else

	override AUTOCLEAN=false
	CONTINUE=true
endif # END IF BUILT CORE CHECK

ifneq ($(wildcard namelist.$(NAMELIST_SUFFIX)), ) # Check for generated namelist file.
	NAMELIST_MESSAGE="A default namelist file (namelist.$(NAMELIST_SUFFIX).defaults) has been generated, but namelist.$(NAMELIST_SUFFIX) has not been modified."
else
	NAMELIST_MESSAGE="A default namelist file (namelist.$(NAMELIST_SUFFIX).defaults) has been generated and copied to namelist.$(NAMELIST_SUFFIX)."
endif


ifeq "$(findstring clean, $(MAKECMDGOALS))" "clean" # CHECK FOR CLEAN TARGET
	override AUTOCLEAN=false
endif # END OF CLEAN TARGET CHECK

VER=$(shell git describe --dirty 2> /dev/null)
#override CPPFLAGS += -DMPAS_GIT_VERSION=$(VER)

ifeq "$(findstring v, $(VER))" "v"
	override CPPFLAGS += -DMPAS_GIT_VERSION=$(VER)
else
	override CPPFLAGS += -DMPAS_GIT_VERSION="unknown"
endif # END OF GIT DESCRIBE VERSION

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

ifeq ($(wildcard src/core_$(CORE)), ) # CHECK FOR EXISTENCE OF CORE DIRECTORY

all: core_error

else

ifeq ($(wildcard src/core_$(CORE)/build_options.mk), ) # Check for build_options.mk
report_builds:
	@echo "CORE=$(CORE)"
endif

ifeq "$(CONTINUE)" "true"
all: mpas_main
else
all: clean_core
endif

endif


mpas_main:
ifeq "$(AUTOCLEAN)" "true"
	$(RM) .mpas_core_*
endif
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
                 CORE="$(CORE)"\
                 AUTOCLEAN="$(AUTOCLEAN)" \
                 GEN_F90="$(GEN_F90)" \
                 NAMELIST_SUFFIX="$(NAMELIST_SUFFIX)" \
                 EXE_NAME="$(EXE_NAME)"

	@echo "$(EXE_NAME)" > .mpas_core_$(CORE)
	if [ -e src/$(EXE_NAME) ]; then mv src/$(EXE_NAME) .; fi
	if [ -e src/inc/namelist.$(NAMELIST_SUFFIX).defaults ]; then mv src/inc/namelist.$(NAMELIST_SUFFIX).defaults .; fi
	if [ ! -e namelist.$(NAMELIST_SUFFIX) ]; then cp namelist.$(NAMELIST_SUFFIX).defaults namelist.$(NAMELIST_SUFFIX); fi
	@echo "*******************************************************************************"
	@echo $(DEBUG_MESSAGE)
	@echo $(PARALLEL_MESSAGE)
	@echo $(PAPI_MESSAGE)
	@echo $(TAU_MESSAGE)
ifeq "$(AUTOCLEAN)" "true"
	@echo $(AUTOCLEAN_MESSAGE)
endif
	@echo $(GEN_F90_MESSAGE)
	@echo $(NAMELIST_MESSAGE)
	@echo "*******************************************************************************"
clean:
	cd src; $(MAKE) clean RM="$(RM)" CORE="$(CORE)"
	$(RM) .mpas_core_*
	$(RM) $(EXE_NAME)
	$(RM) namelist.$(NAMELIST_SUFFIX).defaults
core_error:
	@echo ""
	@echo "*******************************************************************************"
	@echo "     The directory src/core_$(CORE) does not exist."
	@echo "     $(CORE) is not a valid core choice."
	@echo "*******************************************************************************"
	@echo ""
	exit 1
error: errmsg

clean_core:
	@echo ""
	@echo "*******************************************************************************"
	@echo " The MPAS infrastructure is currently built for the $(LAST_CORE) core."
	@echo " Before building the $(CORE) core, please do one of the following."
	@echo ""
	@echo ""
	@echo " To remove the $(LAST_CORE)_model executable and clean the MPAS infrastructure, run:"
	@echo "      make clean CORE=$(LAST_CORE)"
	@echo ""
	@echo " To preserve all executables except $(CORE)_model and clean the MPAS infrastructure, run:"
	@echo "      make clean CORE=$(CORE)"
	@echo ""
	@echo " Alternatively, AUTOCLEAN=true can be appended to the make command to force a clean,"
	@echo " build a new $(CORE)_model executable, and preserve all other executables."
	@echo ""
	@echo "*******************************************************************************"
	@echo ""
	exit 1

else # CORE IF

all: error
clean: errmsg
	exit 1
error: errmsg
	@echo "************ ERROR ************"
	@echo "No CORE specified. Quitting."
	@echo "************ ERROR ************"
	@echo ""
	exit 1

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
	@echo "    DEBUG=true    - builds debug version. Default is optimized version."
	@echo "    USE_PAPI=true - builds version using PAPI for timers. Default is off."
	@echo "    TAU=true      - builds version using TAU hooks for profiling. Default is off."
	@echo "    AUTOCLEAN=true    - forces a clean of infrastructure prior to build new core."
	@echo "    GEN_F90=true  - Generates intermediate .f90 files through CPP, and builds with them."
	@echo ""
	@echo "Ensure that NETCDF, PNETCDF, PIO, and PAPI (if USE_PAPI=true) are environment variables"
	@echo "that point to the absolute paths for the libraries."
	@echo ""
ifdef CORE
	exit 1
endif

