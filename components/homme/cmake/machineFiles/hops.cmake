# copied from weaver
## config that worked for this file:
#
#  1) intel/21.3.0                 7) aue/hdf5/1.14.2-oneapi-2024.1.0-openmpi-4.1.6             13) melcor/gcc/13.1.0
#  2) openmpi-intel/4.1.6-cuda     8) aue/openmpi/4.1.6-oneapi-2024.1.0                         14) melcor/git/2.43.0
#  3) aue/binutils/2.41            9) aue/parallel-netcdf/1.12.3-oneapi-2024.1.0-openmpi-4.1.6  15) melcor/python/3.11.5
#  4) aue/gcc/10.3.0              10) aue/netcdf-c/4.9.2-oneapi-2024.1.0-openmpi-4.1.6          16) totalview/2023.3.10
#  5) aue/cuda/11.8.0-gcc-10.3.0  11) aue/cmake/3.27.7                                          17) melcor/gfortran
#  6) aue/zlib/1.3                12) melcor/cmake/3.27.8                                       18) cudatoolkit/12.4
 

#EXECUTE_PROCESS(COMMAND nf-config --prefix
#  RESULT_VARIABLE NFCONFIG_RESULT
#  OUTPUT_VARIABLE NFCONFIG_OUTPUT
#  ERROR_VARIABLE  NFCONFIG_ERROR
#  OUTPUT_STRIP_TRAILING_WHITESPACE
#)
#SET (NetCDF_Fortran_PATH "${NFCONFIG_OUTPUT}" CACHE STRING "")

#EXECUTE_PROCESS(COMMAND nc-config --prefix
#  RESULT_VARIABLE NCCONFIG_RESULT
#  OUTPUT_VARIABLE NCCONFIG_OUTPUT
#  ERROR_VARIABLE  NCCONFIG_ERROR
#  OUTPUT_STRIP_TRAILING_WHITESPACE
#)
#SET (NetCDF_C_PATH "${NCCONFIG_OUTPUT}" CACHE STRING "")


#SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")
SET(BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")

# TPL settings
set(HDF5_DIR $ENV{HDF5_ROOT} CACHE FILEPATH "")
set(ZLIB_DIR $ENV{ZLIB_ROOT} CACHE FILEPATH "")
set(CURL_ROOT $ENV{CURL_ROOT} CACHE FILEPATH "")
set(CURL_LIBRARY "-L$ENV{CURL_ROOT}/lib -lcurl" CACHE STRING "")

# Flag tweaks
#set(CMAKE_C_FLAGS "-w" CACHE STRING "")
#set(ADD_CXX_FLAGS "-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -Wno-unknown-pragmas -I/projects/ppc64le-pwr9-rhel8/tpls/openmpi/4.1.4/gcc/11.3.0/base/vu2aei6/include" CACHE STRING "")
#set(ADD_Fortran_FLAGS " -I/projects/ppc64le-pwr9-rhel8/tpls/openmpi/4.1.4/gcc/11.3.0/base/vu2aei6/include" CACHE STRING "")
#set(CMAKE_EXE_LINKER_FLAGS "-ldl -lopenblas" CACHE STRING "")
set(OPT_FLAGS "-O3" CACHE STRING "")
set(DEBUG_FLAGS "-ffp-contract=off -g" CACHE STRING "")

# Homme settings
#set(HOMMEXX_VECTOR_SIZE 1 CACHE STRING "")
set(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
set(USE_QUEUING FALSE CACHE BOOL "")
set(USE_TRILINOS OFF CACHE BOOL "")
#set(WITH_PNETCDF FALSE CACHE FILEPATH "")
#set(USE_NUM_PROCS 4 CACHE STRING "Num mpiprocs to use")
#set(HAVE_EXTRAE TRUE CACHE BOOL "")
#set(ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
#set(ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")
#set(HOMMEXX_EXEC_SPACE "CUDA" CACHE STRING "")
#set(HOMMEXX_CUDA_MAX_WARP_PER_TEAM 8 CACHE STRING "")

#SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")
SET(HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")

SET (CPRNC_DIR /projects/ccsm/cprnc/build.toss3 CACHE FILEPATH "")

# Kokkos settings
#set(ENABLE_OPENMP FALSE CACHE BOOL "")
#set(Kokkos_ENABLE_DEBUG FALSE CACHE BOOL "")
#set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION FALSE CACHE BOOL "")
#set(Kokkos_ENABLE_CUDA TRUE CACHE BOOL "")
#set(Kokkos_ENABLE_CUDA_LAMBDA TRUE CACHE BOOL "")
#set(Kokkos_ARCH_VOLTA70 TRUE CACHE BOOL "")
#set(Kokkos_ENABLE_DEPRECATED_CODE FALSE CACHE BOOL "")
#set(Kokkos_ENABLE_EXPLICIT_INSTANTIATION FALSE CACHE BOOL "")

# Compilers
set(CMAKE_C_COMPILER mpicc CACHE STRING "")
#set(CMAKE_CXX_COMPILER mpicxx CACHE STRING "")
set(CMAKE_Fortran_COMPILER mpif90 CACHE STRING "")
######## CHANGE THAT
set(CMAKE_CXX_COMPILER "/home/saruite/E3SM/externals/ekat/extern/kokkos/bin/nvcc_wrapper" CACHE STRING "")
