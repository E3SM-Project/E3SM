# Bit-For-Bit Testing CMake Cache File for Weaver (IBM POWER 9 & NVIDIA TESLA V100s)
# Skyler Ruiter -- 5/30/24
#
#  Module list that worked for this cache file (throw into a script load-weaver)
#
#  1) cmake/3.25.1    4) py-netcdf4/1.5.8   7) openmpi/4.1.4         10) parallel-netcdf/1.12.3  13) curl/7.87.0
#  2) git/2.39.1      5) gcc/11.3.0         8) netcdf-c/4.9.0        11) netlib-lapack/3.10.1
#  3) python/3.10.8   6) cuda/11.8.0        9) netcdf-fortran/4.6.0  12) openblas/0.3.23


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
SET (NetCDF_Fortran_PATH "${NFCONFIG_OUTPUT}" CACHE STRING "") # NEEDED

EXECUTE_PROCESS(COMMAND nc-config --prefix
  RESULT_VARIABLE NCCONFIG_RESULT
  OUTPUT_VARIABLE NCCONFIG_OUTPUT
  ERROR_VARIABLE  NCCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
SET (NetCDF_C_PATH "${NCCONFIG_OUTPUT}" CACHE STRING "") # NEEDED

# TPL settings
#SET (PnetCDF_PATH /projects/ppc64le-pwr9-rhel8/tpls/parallel-netcdf/1.12.3/gcc/11.3.0/openmpi/4.1.4/pp6jcjk CACHE FILEPATH "")
SET (CPRNC_DIR /projects/e3sm/cprnc CACHE FILEPATH "")

# Compilers
SET (CMAKE_C_COMPILER mpicc CACHE STRING "")
SET (CMAKE_Fortran_COMPILER mpif90 CACHE STRING "")
SET (CMAKE_CXX_COMPILER "${CMAKE_CURRENT_SOURCE_DIR}/../../externals/ekat/extern/kokkos/bin/nvcc_wrapper" CACHE STRING "")


############################################
#                                          #
#             BUILD SETTINGS               #
#                                          #
############################################

# Flag settings
SET (CMAKE_C_FLAGS "-w" CACHE STRING "") # disable warnings
SET (ADD_CXX_FLAGS "-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -Wno-unknown-pragmas -I/projects/ppc64le-pwr9-rhel8/tpls/openmpi/4.1.4/gcc/11.3.0/base/vu2aei6/include" CACHE STRING "")
SET (ADD_Fortran_FLAGS " -I/projects/ppc64le-pwr9-rhel8/tpls/openmpi/4.1.4/gcc/11.3.0/base/vu2aei6/include" CACHE STRING "")
SET (CMAKE_EXE_LINKER_FLAGS "-ldl -lopenblas" CACHE STRING "")
SET (OPT_FLAGS "-O0" CACHE STRING "")
SET (DEBUG_FLAGS "-ffp-contract=off -g" CACHE STRING "")

# Homme settings
SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")			# Bit-For-Bit Testing
SET (BUILD_HOMME_WITHOUT_PIOLIBRARY FALSE CACHE BOOL "")	# Buidling without Parallel IO Lib
SET (HOMMEXX_VECTOR_SIZE 1 CACHE STRING "")			# Vector size
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")			# Have homme BLAS
SET (USE_QUEUING FALSE CACHE BOOL "")				# Use queing or not
SET (USE_TRILINOS OFF CACHE BOOL "")				# Use trilinos or not
SET (WITH_PNETCDF FALSE CACHE FILEPATH "")			# Use parallel netcdf
SET (USE_NUM_PROCS 4 CACHE STRING "Num mpiprocs to use")	# Num MPI processes to use
SET (HAVE_EXTRAE TRUE CACHE BOOL "")				# 
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


##########################
#                        #
#      CMake Notes       #
#                        #
##########################
#
# - running make on the individual build gives a few warnings about kokkos
#   things that have changed since newer versions of gcc came out (5 specifically)
# 
# - warnings when running cmake with the cache file: these are about names being
#   being passed to functions being slightly off. Not sure how to remove, see:
#     - https://github.com/zeromq/zproject/pull/1271
#     - https://github.com/zeromq/czmq/issues/2130
#     - https://github.com/gnuradio/gnuradio/issues/3581
#
# - Also error for policy: Policy CMP0075 is not set: Include file check macros honor
#   CMAKE_REQUIRED_LIBRARIES.  Run "cmake --help-policy CMP0075" for policy
#   details.  Use the cmake_policy command to set the policy and suppress this warning.
#   Also not sure how to fix.
