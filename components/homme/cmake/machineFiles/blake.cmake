# Performance Cache File for Blake
# copied from hops
## config that worked for this file:
#
# Currently Loaded Modules:
#   1) gcc/11.3.0     3) openblas/0.3.23   5) openmpi/4.1.5          7) netlib-lapack/3.11.0   9) curl/8.1.2      11) parallel-netcdf/1.12.3
#   2) cmake/3.27.4   4) cuda/11.8.0       6) netcdf-fortran/4.6.0   8) git/2.41.0            10) netcdf-c/4.9.2

# Blake uses Slurm:
# 1) run cmake
# 2) make your builds
# 3) get a node (ex: salloc -N1 --time=00:10:00 --account=fyXXXXXX --partition=H100)
# 	PS- you can run mywcid for the account num and make sure partition is for H100 node
# 4) run mpiexec for your build

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


SET(HOMMEXX_BFB_TESTING FALSE CACHE BOOL "")
SET(BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")

# TPL settings
SET (CPRNC_DIR /projects/e3sm/cprnc CACHE FILEPATH "")

# Flag tweaks
SET (CMAKE_C_FLAGS "-w" CACHE STRING "") # disable warnings
SET (ADD_CXX_FLAGS "-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -Wno-unknown-pragmas -I/projects/ppc64le-pwr9-rhel8/tpls/openmpi/4.1.4/gcc/11.3.0/base/vu2aei6/include" CACHE STRING "")
SET (ADD_Fortran_FLAGS " -I/projects/ppc64le-pwr9-rhel8/tpls/openmpi/4.1.4/gcc/11.3.0/base/vu2aei6/include" CACHE STRING "")
SET (CMAKE_EXE_LINKER_FLAGS "-ldl -lopenblas" CACHE STRING "")
SET (OPT_FLAGS "-O3" CACHE STRING "")
SET (DEBUG_FLAGS "-ffp-contract=off -g" CACHE STRING "")

# Homme settings
set(HOMMEXX_VECTOR_SIZE 1 CACHE STRING "")
set(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
set(USE_QUEUING FALSE CACHE BOOL "")
set(USE_TRILINOS OFF CACHE BOOL "")
set(WITH_PNETCDF FALSE CACHE FILEPATH "")
set(USE_NUM_PROCS 2 CACHE STRING "Num mpiprocs to use")
set(HAVE_EXTRAE TRUE CACHE BOOL "")
set(ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set(ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")
set(HOMMEXX_EXEC_SPACE "CUDA" CACHE STRING "")
set(HOMMEXX_CUDA_MAX_WARP_PER_TEAM 8 CACHE STRING "")
SET(HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")
#SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")
SET(HOMME_ENABLE_COMPOSE TRUE CACHE BOOL "")

# Kokkos settings
set(ENABLE_OPENMP FALSE CACHE BOOL "")
set(Kokkos_ENABLE_DEBUG FALSE CACHE BOOL "")
set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION FALSE CACHE BOOL "")
set(Kokkos_ENABLE_CUDA TRUE CACHE BOOL "")
set(Kokkos_ENABLE_CUDA_LAMBDA TRUE CACHE BOOL "")
set(Kokkos_ARCH_HOPPER90 TRUE CACHE BOOL "") # H100 Architecture Setting
set(Kokkos_ENABLE_DEPRECATED_CODE FALSE CACHE BOOL "")
set(Kokkos_ENABLE_EXPLICIT_INSTANTIATION FALSE CACHE BOOL "")

# Compilers
SET(CMAKE_C_COMPILER mpicc CACHE STRING "")
SET(CMAKE_Fortran_COMPILER mpif90 CACHE STRING "")
SET(CMAKE_CXX_COMPILER "${CMAKE_CURRENT_SOURCE_DIR}/../../externals/ekat/extern/kokkos/bin/nvcc_wrapper" CACHE STRING "")
