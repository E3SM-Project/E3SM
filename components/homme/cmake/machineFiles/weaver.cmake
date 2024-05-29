# CMake Cache File for Weaver (IBM POWER 9 & NVIDIA TESLA V100s)
#
#  Module list that worked for this cache file (throw into a script load-weaver)
#
#  1) binutils/2.30.0                                       8) curl/7.46.0
#  2) gcc/8.3.1                                             9) openblas/0.3.18/gcc/8.3.1
#  3) cuda/11.2.2/gcc/8.3.1                                10) hdf5/1.10.7/gcc/8.3.1/openmpi/4.1.1
#  4) numa/2.0.11                                          11) netcdf-c/4.8.1/gcc/8.3.1/openmpi/4.1.1
#  5) openmpi/4.1.1/gcc/8.3.1/cuda/11.2.2                  12) netcdf-cxx/4.2/gcc/8.3.1/openmpi/4.1.1
#  6) openssl/1.0.2k                                       14) netcdf-fortran/4.5.4/gcc/8.3.1/openmpi/4.1.1
#  7) zlib/1.2.8                                           15) parallel-netcdf/1.12.2/gcc/8.3.1/openmpi/4.1.1
#                                                          16) cmake/3.23.1


############################################
#                                          #
#                ENV PATHS                 #
#                                          #
############################################

# PIO settings
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

# TPL settings
SET (HDF5_DIR $ENV{HDF5_ROOT} CACHE FILEPATH "")
SET (ZLIB_DIR $ENV{ZLIB_ROOT} CACHE FILEPATH "")
SET (CURL_ROOT $ENV{CURL_ROOT} CACHE FILEPATH "")
SET (CURL_LIBRARY "-L$ENV{CURL_ROOT}/lib -lcurl" CACHE STRING "")

SET (CPRNC_DIR /projects/e3sm/cprnc CACHE FILEPATH "")

# Compilers
SET (CMAKE_C_COMPILER mpicc CACHE STRING "")
#SET (CMAKE_CXX_COMPILER mpicxx CACHE STRING "")
SET (CMAKE_Fortran_COMPILER mpif90 CACHE STRING "")
SET (CMAKE_CXX_COMPILER "/home/saruite/E3SM/externals/ekat/extern/kokkos/bin/nvcc_wrapper" CACHE STRING "")


############################################
#                                          #
#             BUILD SETTINGS               #
#                                          #
############################################

# Flag settings
SET (CMAKE_C_FLAGS "-w" CACHE STRING "")
SET (ADD_CXX_FLAGS "-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -Wno-unknown-pragmas -I/projects/ppc64le-pwr9-rhel8/tpls/openmpi/4.1.4/gcc/11.3.0/base/vu2aei6/include" CACHE STRING "")
SET (ADD_Fortran_FLAGS " -I/projects/ppc64le-pwr9-rhel8/tpls/openmpi/4.1.4/gcc/11.3.0/base/vu2aei6/include" CACHE STRING "")
SET (CMAKE_EXE_LINKER_FLAGS "-ldl -lopenblas" CACHE STRING "")
SET (OPT_FLAGS "-O3" CACHE STRING "")
SET (DEBUG_FLAGS "-ffp-contract=off -g" CACHE STRING "")

# Homme settings
#SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")			# Bit-For-Bit Testing
SET (BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")		# Buidling without Parallel IO Lib
SET (HOMMEXX_VECTOR_SIZE 1 CACHE STRING "")			# Vector size
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")			# Have homme BLAS
SET (USE_QUEUING FALSE CACHE BOOL "")				# Use queing or not
SET (USE_TRILINOS OFF CACHE BOOL "")				# Use trilinos or not
#SET (WITH_PNETCDF FALSE CACHE FILEPATH "")			# Use parallel netcdf
SET (USE_NUM_PROCS 4 CACHE STRING "Num mpiprocs to use")	# Num MPI processes to use
SET (HAVE_EXTRAE TRUE CACHE BOOL "")				# idk
SET (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")			# Use column orientation for openmp
SET (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")			# Use horizontal orientation for openmp
SET (HOMMEXX_EXEC_SPACE "CUDA" CACHE STRING "")			# Define the execution space for homme
SET (HOMMEXX_CUDA_MAX_WARP_PER_TEAM 8 CACHE STRING "")		# Cuda warps per team setting
#SET (BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")		# Build preqx with kokkos
SET (BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")		# Build theta with kokkos
SET (HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")			# Enabling compose

# Kokkos settings
SET (ENABLE_OPENMP FALSE CACHE BOOL "")				# 
SET (Kokkos_ENABLE_DEBUG FALSE CACHE BOOL "")			# Debug flag
SET (Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION FALSE CACHE BOOL "")#
SET (Kokkos_ENABLE_CUDA TRUE CACHE BOOL "")			# Enable Cuda
SET (Kokkos_ENABLE_CUDA_LAMBDA TRUE CACHE BOOL "")		# Enable cuda lambdas
SET (Kokkos_ARCH_VOLTA70 TRUE CACHE BOOL "")			# Informing kokkos on the hardware
SET (Kokkos_ENABLE_DEPRECATED_CODE FALSE CACHE BOOL "")		# Enable deprecated code (Sometimes needed?)
SET (Kokkos_ENABLE_EXPLICIT_INSTANTIATION FALSE CACHE BOOL "")	# Enable explicit instantiation
