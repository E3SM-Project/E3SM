## config that worked for this file:
## export PATH=/ascldap/users/projects/e3sm/scream/libs/netcdf-fortran/install/weaver/bin:$PATH
#
#  1) binutils/2.30.0                                       8) curl/7.46.0
#  2) gcc/7.2.0                                             9) openblas/0.2.20/gcc/7.2.0
#  3) cuda/10.1.105                                        10) hdf5/1.10.5/openmpi/4.0.1/gcc/7.2.0/cuda/10.1.105
#  4) numa/2.0.11                                          11) pnetcdf/1.9.0/openmpi/4.0.1/gcc/7.2.0/cuda/10.1.105
#  5) openmpi/4.0.1/gcc/7.2.0/cuda/10.1.105                12) netcdf/4.6.1/openmpi/4.0.1/gcc/7.2.0/cuda/10.1.105
#  6) openssl/1.0.2k                                       13) cmake/3.19.3
#  7) zlib/1.2.8



EXECUTE_PROCESS(COMMAND nf-config --prefix
  RESULT_VARIABLE NFCONFIG_RESULT
  OUTPUT_VARIABLE NFCONFIG_OUTPUT
  ERROR_VARIABLE  NFCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
SET (NetCDF_Fortran_PATH "${NFCONFIG_OUTPUT}" CACHE STRING "")

EXECUTE_PROCESS(COMMAND nc-config --prefix
  RESULT_VARIABLE NCCONFIG_RESULT
  OUTPUT_VARIABLE NCCONFIG_OUTPUT
  ERROR_VARIABLE  NCCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
SET (NetCDF_C_PATH "${NCCONFIG_OUTPUT}" CACHE STRING "")


SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

# TPL settings
set(HDF5_DIR $ENV{HDF5_ROOT} CACHE FILEPATH "")
set(ZLIB_DIR $ENV{ZLIB_ROOT} CACHE FILEPATH "")
set(CURL_ROOT $ENV{CURL_ROOT} CACHE FILEPATH "")
set(CURL_LIBRARY -L$ENV{CURL_ROOT}/lib -lcurl CACHE LIST "")

# Flag tweaks
set(CMAKE_C_FLAGS "-w" CACHE STRING "")
set(ADD_CXX_FLAGS "-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -Wno-unknown-pragmas --fmad=false -O0" CACHE STRING "")
set(ADD_Fortran_FLAGS " -ffp-contract=off -O0" CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS "-ldl" CACHE STRING "")
set(OPT_FLAGS "-O0" CACHE STRING "")
set(DEBUG_FLAGS "-ffp-contract=off -g"CACHE STRING "")

# Homme settings
set(HOMMEXX_VECTOR_SIZE 1 CACHE STRING "")
set(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
set(USE_QUEUING FALSE CACHE BOOL "")
set(USE_TRILINOS OFF CACHE BOOL "")
#set(WITH_PNETCDF FALSE CACHE FILEPATH "")
set(USE_NUM_PROCS 4 CACHE STRING "Num mpiprocs to use")
set(HAVE_EXTRAE TRUE CACHE BOOL "")
set(ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set(ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")
set(HOMMEXX_EXEC_SPACE "CUDA" CACHE STRING "")


SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")
#sET(HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")

set(HOMMEXX_BFB_TESTING TRUE CACHE STRING ON)

# Kokkos settings
set(ENABLE_OPENMP FALSE CACHE BOOL "")
set(Kokkos_ENABLE_DEBUG FALSE CACHE BOOL "")
set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION FALSE CACHE BOOL "")
set(Kokkos_ENABLE_CUDA TRUE CACHE BOOL "")
set(Kokkos_ENABLE_CUDA_LAMBDA TRUE CACHE BOOL "")
set(Kokkos_ARCH_VOLTA70 TRUE CACHE BOOL "")
set(Kokkos_ENABLE_DEPRECATED_CODE FALSE CACHE BOOL "")
set(Kokkos_ENABLE_EXPLICIT_INSTANTIATION FALSE CACHE BOOL "")

# Compilers
set(CMAKE_C_COMPILER mpicc CACHE STRING "")
set(CMAKE_Fortran_COMPILER mpif90 CACHE STRING "")
######## CHANGE THAT
set(CMAKE_CXX_COMPILER "/home/onguba/acme-master/externals/kokkos/bin/nvcc_wrapper" CACHE STRING "")
