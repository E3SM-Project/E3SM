# CMake initial cache file for the ghci gnu image, for BFB (cxx vs f90) testing
#
# NOTE: these must be set BEFORE including the non-bfb file (the first value set in the cache wins)
SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

# Turn off optimizations and fma contractions, so that Fortran and C++ results can match bit-for-bit
SET (ADD_Fortran_FLAGS "-fbacktrace -O0 -ffp-contract=off" CACHE STRING "")
SET (ADD_C_FLAGS       "-O0 -ffp-contract=off" CACHE STRING "")
SET (ADD_CXX_FLAGS     "-O0 -ffp-contract=off" CACHE STRING "")

include (${CMAKE_CURRENT_LIST_DIR}/ghci-gnu.cmake)
