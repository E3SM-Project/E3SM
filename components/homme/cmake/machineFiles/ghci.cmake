# CMake initial cache file for our GitHub CI (ghci) self-hosted runners
#
# Common settings for all ghci images. Compiler-specific settings
# are in ghci-gnu.cmake, ghci-intel.cmake, and ghci-cuda.cmake
# (the latter are picked up by test_homme.py via the --machine-file flag)
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER       mpicc  CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER     mpicxx CACHE FILEPATH "")

SET (WITH_PNETCDF FALSE CACHE BOOL "")
SET (NetCDF_Fortran_PATH "$ENV{NETCDF_FORTRAN_ROOT}" CACHE STRING "" FORCE)
SET (NetCDF_C_PATH       "$ENV{NETCDF_C_ROOT}"       CACHE STRING "" FORCE)

# NOTE: we do not set CPRNC_DIR, since there is no prebuilt cprnc in our images. Homme will compile it
# from the CIME sources instead (it is built along with the test executables).

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (USE_MPIEXEC "mpirun" CACHE STRING "")

SET (HOMMEXX_BFB_TESTING      FALSE   CACHE BOOL "")
SET (BUILD_HOMME_PREQX_KOKKOS TRUE    CACHE BOOL "")
SET (BUILD_HOMME_THETA_KOKKOS TRUE    CACHE BOOL "")
SET (HOMME_TESTING_PROFILE    "short" CACHE STRING "")
SET (HOMME_TESTING_TIMELIMIT  1800    CACHE STRING "")

SET (BUILD_HOMME_WITHOUT_PIOLIBRARY FALSE CACHE BOOL "")
SET (USE_TRILINOS OFF CACHE BOOL "")

# Let's catch usage of code deprecated in Kokkos 4
SET (Kokkos_ENABLE_DEPRECATED_CODE_4 OFF CACHE BOOL "")
# NOTE: we must keep the code deprecated in Kokkos 5 enabled (which is also the Kokkos default): without it,
# Kokkos aborts if the vector length of a TeamPolicy is not a power of 2, and Homme's CPU code uses the number
# of level packs (e.g., ceil(20/8)=3 or ceil(72/8)=9) as vector length in some TeamPolicy (BoundaryExchange.cpp).
SET (Kokkos_ENABLE_DEPRECATED_CODE_5 ON CACHE BOOL "")
