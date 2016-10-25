MODEL_FORMULATION = 


dummy:
	( $(MAKE) error )

xlf:
	( $(MAKE) all \
	"FC_PARALLEL = mpxlf90" \
	"CC_PARALLEL = mpcc" \
	"CXX_PARALLEL = mpixlcxx" \
	"FC_SERIAL = xlf90" \
	"CC_SERIAL = xlc" \
	"CXX_SERIAL = xlcxx" \
	"FFLAGS_PROMOTION = -qrealsize=8" \
	"FFLAGS_OPT = -O3" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -O0 -g -C" \
	"CFLAGS_DEBUG = -O0 -g" \
	"CXXFLAGS_DEBUG = -O0 -g" \
	"LDFLAGS_DEBUG = -O0 -g" \
	"FFLAGS_OMP = -qsmp=omp" \
	"CFLAGS_OMP = -qsmp=omp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI" )
 
ftn:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"CXX_PARALLEL = CC" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"CXX_SERIAL = CC" \
	"FFLAGS_PROMOTION = -r8" \
	"FFLAGS_OPT = -i4 -gopt -O2 -Mvect=nosse -Kieee -convert big_endian" \
	"CFLAGS_OPT = -fast" \
	"CXXFLAGS_OPT = -fast" \
	"LDFLAGS_OPT = " \
	"FFLAGS_OMP = -mp" \
	"CFLAGS_OMP = -mp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

titan-cray:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = gcc" \
	"FFLAGS_PROMOTION = -default64" \
	"FFLAGS_OPT = -s integer32 -O3 -f free -N 255 -em -ef" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_OMP = " \
	"CFLAGS_OMP = " \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

pgi:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"CXX_PARALLEL = mpicxx" \
	"FC_SERIAL = pgf90" \
	"CC_SERIAL = pgcc" \
	"CXX_SERIAL = pgc++" \
	"FFLAGS_PROMOTION = -r8" \
	"FFLAGS_OPT = -O3 -byteswapio -Mfree" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -O0 -g -Mbounds -Mchkptr -byteswapio -Mfree -Ktrap=divz,fp,inv,ovf -traceback" \
	"CFLAGS_DEBUG = -O0 -g -traceback" \
	"CXXFLAGS_DEBUG = -O0 -g -traceback" \
	"LDFLAGS_DEBUG = -O0 -g -Mbounds -Mchkptr -Ktrap=divz,fp,inv,ovf -traceback" \
	"FFLAGS_OMP = -mp" \
	"CFLAGS_OMP = -mp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

pgi-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"CXX_PARALLEL = CC" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"CXX_SERIAL = CC" \
	"FFLAGS_PROMOTION = -r8" \
	"FFLAGS_OPT = -O3 -byteswapio -Mfree" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_OMP = -mp" \
	"CFLAGS_OMP = -mp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

pgi-llnl:
	( $(MAKE) all \
	"FC_PARALLEL = mpipgf90" \
	"CC_PARALLEL = pgcc" \
	"CXX_PARALLEL = mpipgcxx" \
	"FC_SERIAL = pgf90" \
	"CC_SERIAL = pgcc" \
	"CXX_SERIAL = pgc++" \
	"FFLAGS_PROMOTION = -r8" \
	"FFLAGS_OPT = -i4 -g -O2 -byteswapio" \
	"CFLAGS_OPT = -fast" \
	"CXXFLAGS_OPT = -fast" \
	"LDFLAGS_OPT = " \
	"FFLAGS_OMP = -mp" \
	"CFLAGS_OMP = -mp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

ifort:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"CXX_PARALLEL = mpicxx" \
	"FC_SERIAL = ifort" \
	"CC_SERIAL = icc" \
	"CXX_SERIAL = icpc" \
	"FFLAGS_PROMOTION = -real-size 64" \
	"FFLAGS_OPT = -O3 -convert big_endian -FR" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -g -convert big_endian -FR -CU -CB -check all -fpe0 -traceback" \
	"CFLAGS_DEBUG = -g -fpe0 -traceback" \
	"CXXFLAGS_DEBUG = -g -fpe0 -traceback" \
	"LDFLAGS_DEBUG = -g -fpe0 -traceback" \
	"FFLAGS_OMP = -qopenmp" \
	"CFLAGS_OMP = -qopenmp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

ifort-mic:
	( $(MAKE) all \
	"FC_PARALLEL = mpiifort" \
	"CC_PARALLEL = mpiicc" \
	"CXX_PARALLEL = mpiicpc" \
	"FC_SERIAL = ifort" \
	"CC_SERIAL = icc" \
	"CXX_SERIAL = icpc" \
	"FFLAGS_PROMOTION = -real-size 64" \
	"FFLAGS_OPT = -O2 -mmic -convert big_endian -FR" \
	"CFLAGS_OPT = -O2 -mmic" \
	"CXXFLAGS_OPT = -O2 -mmic" \
	"LDFLAGS_OPT = -O2 -mmic" \
	"FFLAGS_DEBUG = -O0 -mmic -g -convert big_endian -FR -CU -CB -check all -fp-model strict -traceback" \
	"CFLAGS_DEBUG = -O0 -mmic -g -fp-model strict -traceback" \
	"CXXFLAGS_DEBUG = -O0 -mmic -g -fp-model strict -traceback" \
	"LDFLAGS_DEBUG = -O0 -mmic -g -fp-model strict -traceback" \
	"FFLAGS_OMP = -qopenmp" \
	"CFLAGS_OMP = -qopenmp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

ifort-scorep:
	( $(MAKE) all \
	"FC_PARALLEL = scorep --compiler mpif90" \
	"CC_PARALLEL = scorep --compiler mpicc" \
	"CXX_PARALLEL = scorep --compiler mpicxx" \
	"FC_SERIAL = ifort" \
	"CC_SERIAL = icc" \
	"CXX_SERIAL = icpc" \
	"FFLAGS_PROMOTION = -real-size 64" \
	"FFLAGS_OPT = -O3 -convert big_endian -FR" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -g -convert big_endian -FR -CU -CB -check all -fpe0 -traceback" \
	"CFLAGS_DEBUG = -g -fpe0 -traceback" \
	"CXXFLAGS_DEBUG = -g -fpe0 -traceback" \
	"LDFLAGS_DEBUG = -g -fpe0 -traceback" \
	"FFLAGS_OMP = -qopenmp" \
	"CFLAGS_OMP = -qopenmp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

ifort-gcc:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"CXX_PARALLEL = mpicxx" \
	"FC_SERIAL = ifort" \
	"CC_SERIAL = gcc" \
	"CXX_SERIAL = g++" \
	"FFLAGS_PROMOTION = -real-size 64" \
	"FFLAGS_OPT = -O3 -convert big_endian -FR" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -g -convert big_endian -FR -CU -CB -check all -fpe0 -traceback" \
	"CFLAGS_DEBUG = -g" \
	"CXXFLAGS_DEBUG = -g" \
	"LDFLAGS_DEBUG = -g -fpe0 -traceback" \
	"FFLAGS_OMP = -qopenmp" \
	"CFLAGS_OMP = -fopenmp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

gfortran:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"CXX_PARALLEL = mpicxx" \
	"FC_SERIAL = gfortran" \
	"CC_SERIAL = gcc" \
	"CXX_SERIAL = g++" \
	"FFLAGS_PROMOTION = -fdefault-real-8 -fdefault-double-8" \
	"FFLAGS_OPT = -O3 -m64 -ffree-line-length-none -fconvert=big-endian -ffree-form" \
	"CFLAGS_OPT = -O3 -m64" \
	"CXXFLAGS_OPT = -O3 -m64" \
	"LDFLAGS_OPT = -O3 -m64" \
	"FFLAGS_DEBUG = -g -m64 -ffree-line-length-none -fconvert=big-endian -ffree-form -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow" \
	"CFLAGS_DEBUG = -g -m64" \
	"CXXFLAGS_DEBUG = -O3 -m64" \
	"LDFLAGS_DEBUG = -g -m64" \
	"FFLAGS_OMP = -fopenmp" \
	"CFLAGS_OMP = -fopenmp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

gfortran-clang:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc -cc=clang" \
	"CXX_PARALLEL = mpicxx -cxx=clang++" \
	"FC_SERIAL = gfortran" \
	"CC_SERIAL = clang" \
	"CXX_SERIAL = clang++" \
	"FFLAGS_PROMOTION = -fdefault-real-8 -fdefault-double-8" \
	"FFLAGS_OPT = -O3 -m64 -ffree-line-length-none -fconvert=big-endian -ffree-form" \
	"CFLAGS_OPT = -O3 -m64" \
	"CXXFLAGS_OPT = -O3 -m64" \
	"LDFLAGS_OPT = -O3 -m64" \
	"FFLAGS_DEBUG = -g -m64 -ffree-line-length-none -fconvert=big-endian -ffree-form -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow" \
	"CFLAGS_DEBUG = -g -m64" \
	"CXXFLAGS_DEBUG = -O3 -m64" \
	"LDFLAGS_DEBUG = -g -m64" \
	"FFLAGS_OMP = -fopenmp" \
	"CFLAGS_OMP = -fopenmp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

g95:
	( $(MAKE) all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"CXX_PARALLEL = mpicxx" \
	"FC_SERIAL = g95" \
	"CC_SERIAL = gcc" \
	"CXX_SERIAL = g++" \
	"FFLAGS_PROMOTION = -r8" \
	"FFLAGS_OPT = -O3 -ffree-line-length-huge -fendian=big" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_OMP = -fopenmp" \
	"CFLAGS_OMP = -fopenmp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

pathscale-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"CXX_PARALLEL = CC" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"CXX_SERIAL = CC" \
	"FFLAGS_PROMOTION = -r8" \
	"FFLAGS_OPT = -O3 -freeform -extend-source" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_OMP = -mp" \
	"CFLAGS_OMP = -mp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

cray-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"CXX_PARALLEL = CC" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"CXX_SERIAL = CC" \
	"FFLAGS_PROMOTION = -default64" \
	"FFLAGS_OPT = -O3 -f free" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_OMP = " \
	"CFLAGS_OMP = " \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

gnu-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"CXX_PARALLEL = CC" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"CXX_SERIAL = CC" \
	"FFLAGS_PROMOTION = -fdefault-real-8 -fdefault-double-8" \
	"FFLAGS_OPT = -O3 -m64 -ffree-line-length-none -fconvert=big-endian -ffree-form" \
	"CFLAGS_OPT = -O3 -m64" \
	"CXXFLAGS_OPT = -O3 -m64" \
	"LDFLAGS_OPT = -O3 -m64" \
	"FFLAGS_DEBUG = -g -m64 -ffree-line-length-none -fconvert=big-endian -ffree-form" \
	"CFLAGS_DEBUG = -g -m64" \
	"CXXFLAGS_DEBUG = -g -m64" \
	"LDFLAGS_DEBUG = -g -m64" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE -D_MPI $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

intel-nersc:
	( $(MAKE) all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"CXX_PARALLEL = CC" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"CXX_SERIAL = CC" \
	"FFLAGS_PROMOTION = -real-size 64" \
	"FFLAGS_OPT = -O3 -convert big_endian -FR" \
	"CFLAGS_OPT = -O3" \
	"CXXFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_OMP = -qopenmp" \
	"CFLAGS_OMP = -qopenmp" \
	"FFLAGS_DEBUG = -real-size 64 -g -convert big_endian -FR -CU -CB -check all -gen-interfaces -warn interfaces -traceback" \
	"CFLAGS_DEBUG = -g -traceback" \
	"CXXFLAGS_DEBUG = -g -traceback" \
	"LDFLAGS_DEBUG = -g -traceback" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE" )

bluegene:
	( $(MAKE) all \
	"FC_PARALLEL = mpixlf95_r" \
	"CC_PARALLEL = mpixlc_r" \
	"CXX_PARALLEL = mpixlcxx_r" \
	"FC_SERIAL = bgxlf95_r" \
	"CC_SERIAL = bgxlc_r" \
	"CXX_SERIAL = bgxlc++_r" \
	"FFLAGS_PROMOTION = -qrealsize=8" \
	"FFLAGS_OPT = -O2 -g" \
	"CFLAGS_OPT = -O2 -g" \
	"CXXFLAGS_OPT = -O2 -g" \
	"LDFLAGS_OPT = -O2 -g" \
	"FFLAGS_DEBUG = -O0 -g -C -qinitalloc -qinitauto" \
	"CFLAGS_DEBUG = -O0 -g" \
	"CXXFLAGS_DEBUG = -O0 -g" \
	"LDFLAGS_DEBUG = -O0 -g" \
	"FFLAGS_OMP = -qsmp=omp" \
	"CFLAGS_OMP = -qsmp=omp" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"USE_PAPI = $(USE_PAPI)" \
	"OPENMP = $(OPENMP)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -D_MPI" )

CPPINCLUDES = 
FCINCLUDES = 
LIBS = 
ifneq ($(wildcard $(PIO)/lib), ) # Check for newer PIO version
ifeq "$(USE_PIO2)" "true"
	CPPINCLUDES = -DUSE_PIO2 -I$(PIO)/include
	FCINCLUDES = -DUSE_PIO2 -I$(PIO)/include
	LIBS = -L$(PIO)/lib -lpiof -lpioc
else
	CPPINCLUDES = -I$(PIO)/include
	FCINCLUDES = -I$(PIO)/include
	LIBS = -L$(PIO)/lib -lpio
endif
else
ifeq "$(USE_PIO2)" "true"
	CPPINCLUDES = -DUSE_PIO2 -I$(PIO)/include
	FCINCLUDES = -DUSE_PIO2 -I$(PIO)/include
	LIBS = -L$(PIO) -lpiof -lpioc
else
	CPPINCLUDES = -I$(PIO)
	FCINCLUDES = -I$(PIO)
	LIBS = -L$(PIO) -lpio
endif
endif

ifneq "$(PNETCDF)" ""
	CPPINCLUDES += -I$(PNETCDF)/include
	FCINCLUDES += -I$(PNETCDF)/include
	LIBS += -L$(PNETCDF)/lib -lpnetcdf
endif

ifneq "$(NETCDF)" ""
	CPPINCLUDES += -I$(NETCDF)/include
	FCINCLUDES += -I$(NETCDF)/include
	LIBS += -L$(NETCDF)/lib
	NCLIB = -lnetcdf
	NCLIBF = -lnetcdff
	ifneq ($(wildcard $(NETCDF)/lib/libnetcdff.*), ) # CHECK FOR NETCDF4
		LIBS += $(NCLIBF)
	endif # CHECK FOR NETCDF4
	ifneq "$(NETCDFF)" ""
		FCINCLUDES += -I$(NETCDFF)/include
		LIBS += -L$(NETCDFF)/lib
		LIBS += $(NCLIBF)
	endif
	LIBS += $(NCLIB)
endif

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
	CXXFLAGS=$(CXXFLAGS_OPT)
	LDFLAGS=$(LDFLAGS_OPT)
	DEBUG_MESSAGE="Debug flags are not defined for this compile group. Defaulting to Optimized flags"
else # FFLAGS_DEBUG IF
	FFLAGS=$(FFLAGS_DEBUG)
	CFLAGS=$(CFLAGS_DEBUG)
	CXXFLAGS=$(CXXFLAGS_DEBUG)
	LDFLAGS=$(LDFLAGS_DEBUG)
	override CPPFLAGS += -DMPAS_DEBUG
	DEBUG_MESSAGE="Debugging is on."
endif # FFLAGS_DEBUG IF

else # DEBUG IF
	FFLAGS=$(FFLAGS_OPT)
	CFLAGS=$(CFLAGS_OPT)
	CXXFLAGS=$(CXXFLAGS_OPT)
	LDFLAGS=$(LDFLAGS_OPT)
	DEBUG_MESSAGE="Debugging is off."
endif # DEBUG IF

FC=$(FC_PARALLEL)
CC=$(CC_PARALLEL)
CXX=$(CXX_PARALLEL)
SFC=$(FC_SERIAL)
SCC=$(CC_SERIAL)
PARALLEL_MESSAGE="Parallel version is on."

ifeq "$(OPENMP)" "true"
	FFLAGS += $(FFLAGS_OMP)
	CFLAGS += $(CFLAGS_OMP)
	CXXFLAGS += $(CFLAGS_OMP)
	override CPPFLAGS += "-DMPAS_OPENMP"
	LDFLAGS += $(FFLAGS_OMP)
endif #OPENMP IF

ifeq "$(PRECISION)" "single"
	FFLAGS += "-DSINGLE_PRECISION"
	CFLAGS += "-DSINGLE_PRECISION"
	CXXFLAGS += "-DSINGLE_PRECISION"
	override CPPFLAGS += "-DSINGLE_PRECISION"
	PRECISION_MESSAGE="MPAS was built with default single-precision reals."
else
	FFLAGS += $(FFLAGS_PROMOTION)
	PRECISION_MESSAGE="MPAS was built with default double-precision reals."
endif #PRECISION IF

ifeq "$(USE_PAPI)" "true"
	CPPINCLUDES += -I$(PAPI)/include -D_PAPI
	FCINCLUDES += -I$(PAPI)/include
	LIBS += -L$(PAPI)/lib -lpapi
	PAPI_MESSAGE="Papi libraries are on."
else # USE_PAPI IF
	PAPI_MESSAGE="Papi libraries are off."
endif # USE_PAPI IF

ifeq "$(USE_PIO2)" "true"
	PIO_MESSAGE="Using the PIO 2 library."
else # USE_PIO2 IF
	PIO_MESSAGE="Using the PIO 1.x library."
endif # USE_PIO2 IF

ifdef TIMER_LIB
ifeq "$(TIMER_LIB)" "tau"
	override TAU=true
	TIMER_MESSAGE="TAU is being used for the timer interface"
endif

ifeq "$(TIMER_LIB)" "gptl"
	override CPPFLAGS += -DMPAS_GPTL_TIMERS
	override FCINCLUDES += -I${GPTL}/include
	override LIBS += -L${GPTL}/lib -lgptl
	TIMER_MESSAGE="GPTL is being used for the timer interface"
endif

ifeq "$(TIMER_LIB)" ""
	override CPPFLAGS += -DMPAS_NATIVE_TIMERS
	TIMER_MESSAGE="The native timer interface is being used"
endif

else # else ifdef $(TIMER_LIB)

	override CPPFLAGS += -DMPAS_NATIVE_TIMERS
	TIMER_MESSAGE="The native timer interface is being used"

endif # endif ifdef $(TIMER_LIB)

ifeq "$(TAU)" "true"
	LINKER=tau_f90.sh
	CPPINCLUDES += -DMPAS_TAU -DMPAS_TAU_TIMERS
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

ifeq "$(OPENMP)" "true"
	OPENMP_MESSAGE="MPAS was built with OpenMP enabled."
else
	OPENMP_MESSAGE="MPAS was built without OpenMP support."
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

ifneq ($(wildcard streams.$(NAMELIST_SUFFIX)), ) # Check for generated streams file.
	STREAM_MESSAGE="A default streams file (streams.$(NAMELIST_SUFFIX).defaults) has been generated, but streams.$(NAMELIST_SUFFIX) has not been modified."
else
	STREAM_MESSAGE="A default streams file (streams.$(NAMELIST_SUFFIX).defaults) has been generated and copied to streams.$(NAMELIST_SUFFIX)."
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
	override LIBS += $(MPAS_EXTERNAL_LIBS)
endif
ifdef MPAS_EXTERNAL_INCLUDES
	override CPPINCLUDES += $(MPAS_EXTERNAL_INCLUDES)
	override FCINCLUDES += $(MPAS_EXTERNAL_INCLUDES)
endif
ifdef MPAS_EXTERNAL_CPPFLAGS
	override CPPFLAGS += $(MPAS_EXTERNAL_CPPFLAGS)
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


compiler_test:
ifeq "$(OPENMP)" "true"
	@echo "Testing compiler for OpenMP support"
	@echo "int main() { return 0; }" > conftest.c; $(SCC) $(CFLAGS) -o conftest.out conftest.c || \
		(echo "$(SCC) does not support OpenMP - see INSTALL in top-level directory for more information"; rm -fr conftest.*; exit 1)
	@echo "int main() { return 0; }" > conftest.c; $(CC) $(CFLAGS) -o conftest.out conftest.c || \
		(echo "$(CC) does not support OpenMP - see INSTALL in top-level directory for more information"; rm -fr conftest.*; exit 1)
	@echo "int main() { return 0; }" > conftest.cpp; $(CXX) $(CFLAGS) -o conftest.out conftest.cpp || \
		(echo "$(CXX) does not support OpenMP - see INSTALL in top-level directory for more information"; rm -fr conftest.*; exit 1)
	@echo "program test; stop 0; end program" > conftest.f90; $(SFC) $(FFLAGS) -o conftest.out conftest.f90 || \
		(echo "$(SFC) does not support OpenMP - see INSTALL in top-level directory for more information"; rm -fr conftest.*; exit 1)
	@echo "program test; stop 0; end program" > conftest.f90; $(FC) $(FFLAGS) -o conftest.out conftest.f90 || \
		(echo "$(FC) does not support OpenMP - see INSTALL in top-level directory for more information"; rm -fr conftest.*; exit 1)
	@rm -fr conftest.*
endif


mpas_main: compiler_test
ifeq "$(AUTOCLEAN)" "true"
	$(RM) .mpas_core_*
endif
	cd src; $(MAKE) FC="$(FC)" \
                 CC="$(CC)" \
                 CXX="$(CXX)" \
                 SFC="$(SFC)" \
                 SCC="$(SCC)" \
                 LINKER="$(LINKER)" \
                 CFLAGS="$(CFLAGS)" \
                 CXXFLAGS="$(CXXFLAGS)" \
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
	( cd src/core_$(CORE); $(MAKE) ROOT_DIR="$(PWD)" post_build )
	@echo "*******************************************************************************"
	@echo $(PRECISION_MESSAGE)
	@echo $(DEBUG_MESSAGE)
	@echo $(PARALLEL_MESSAGE)
	@echo $(PAPI_MESSAGE)
	@echo $(TAU_MESSAGE)
	@echo $(OPENMP_MESSAGE)
ifeq "$(AUTOCLEAN)" "true"
	@echo $(AUTOCLEAN_MESSAGE)
endif
	@echo $(GEN_F90_MESSAGE)
	@echo $(TIMER_MESSAGE)
	@echo $(PIO_MESSAGE)
	@echo "*******************************************************************************"
clean:
	cd src; $(MAKE) clean RM="$(RM)" CORE="$(CORE)"
	$(RM) .mpas_core_*
	$(RM) $(EXE_NAME)
	$(RM) namelist.$(NAMELIST_SUFFIX).defaults
	$(RM) streams.$(NAMELIST_SUFFIX).defaults
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
	@echo "    TIMER_LIB=opt - Selects the timer library interface to be used for profiling the model. Options are:"
	@echo "                    TIMER_LIB=native - Uses native built-in timers in MPAS"
	@echo "                    TIMER_LIB=gptl - Uses gptl for the timer interface instead of the native interface"
	@echo "                    TIMER_LIB=tau - Uses TAU for the timer interface instead of the native interface"
	@echo "    OPENMP=true   - builds and links with OpenMP flags. Default is to not use OpenMP."
	@echo "    USE_PIO2=true - links with the PIO 2 library. Default is to use the PIO 1.x library."
	@echo "    PRECISION=single - builds with default single-precision real kind. Default is to use double-precision."
	@echo ""
	@echo "Ensure that NETCDF, PNETCDF, PIO, and PAPI (if USE_PAPI=true) are environment variables"
	@echo "that point to the absolute paths for the libraries."
	@echo ""
ifdef CORE
	exit 1
endif

