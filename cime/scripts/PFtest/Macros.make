#
# COMPILER=gnu
# OS=LINUX
# MACH=cades
#
# Makefile Macros 
CPPDEFS+= -DFORTRANUNDERSCORE -DNO_R16 -DCPRGNU  -DFORTRANUNDERSCORE -DNO_R16

BLAS_LIBDIR:=/software/tools/compilers/intel_2017/mkl/lib/intel64

CFLAGS:= -mcmodel=medium 

CXX_LINKER:=FORTRAN

FC_AUTO_R8:= -fdefault-real-8 

FFLAGS:= -O -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -fno-range-check

FFLAGS_NOOPT:= -O0 

FIXEDFLAGS:=  -ffixed-form 

FREEFLAGS:= -ffree-form 

HAS_F2008_CONTIGUOUS:=FALSE

HDF5_PATH:=/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/hdf5-parallel/1.8.17/centos7.2_gnu5.3.0

LAPACK_LIBDIR:=/software/tools/compilers/intel_2017/mkl/lib/intel64

MPICC:=mpicc

MPICXX:=mpic++

MPIFC:=mpif90

NETCDF_PATH:=/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/netcdf-hdf5parallel/4.3.3.1/centos7.2_gnu5.3.0

SCC:=gcc

SCXX:=gcpp

SFC:=gfortran

SUPPORTS_CXX:=TRUE

ifeq ($(DEBUG), FALSE) 
   FFLAGS +=  -O 
   CFLAGS +=  -O 
endif

ifeq ($(DEBUG), TRUE) 
   FFLAGS +=  -g -Wall 
   CFLAGS +=  -g -Wall -Og -fbacktrace -fcheck=bounds -ffpe-trap=invalid,zero,overflow
endif

ifeq ($(compile_threaded), true) 
   FFLAGS +=  -fopenmp 
   LDFLAGS +=  -fopenmp 
   CFLAGS +=  -fopenmp 
endif

ifeq ($(MODEL), cism) 
   CMAKE_OPTS +=  -D CISM_GNU=ON 
endif

ifeq ($(MODEL), clm) 
  ifeq ($(CLM_PFLOTRAN_COLMODE), TRUE) 
    ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
       CPPDEFS +=  -DCOLUMN_MODE 
    endif

  endif

  ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
     FFLAGS +=  -I$(CLM_PFLOTRAN_SOURCE_DIR) 
     CPPDEFS +=  -DCLM_PFLOTRAN 
  endif

endif

ifeq ($(MODEL), csm_share) 
   CFLAGS +=  -std=c99 
endif

ifeq ($(MODEL), driver) 
   LDFLAGS +=  -L$(NETCDF_PATH)/lib -Wl,-rpath=$(NETCDF_PATH)/lib -lnetcdff -lnetcdf \
              -L$(HDF5_PATH)/lib -Wl,-rpath=$(HDF5_PATH)/lib -lhdf5_hl -lhdf5 \
              -L$(LAPACK_LIBDIR) -Wl,-rpath=$(LAPACK_LIBDIR) \
              -L$(BLAS_LIBDIR) -Wl,-rpath=$(BLAS_LIBDIR) \
              -Wl,--start-group -lmkl_blas95_lp64 -lmkl_lapack95_lp64 \
              -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_intel_lp64 \
              -lmkl_core -lmkl_gnu_thread -lmkl_blacs_openmpi_lp64 \
              -Wl,--end-group -lgomp -lrt
      
  ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
     LDFLAGS +=  -L$(CLM_PFLOTRAN_SOURCE_DIR) -lpflotran $(PETSC_LIB) 
  endif

endif

ifeq ($(MODEL), pop) 
   CPPDEFS +=  -D_USE_FLOW_CONTROL 
endif

