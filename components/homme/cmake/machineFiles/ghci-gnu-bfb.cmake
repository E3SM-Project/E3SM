# CMake initial cache file for the ghci gnu image, for BFB (cxx vs f90) testing
#
# NOTE: these must be set BEFORE including the non-bfb file (the first value set in the cache wins)
SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

# Low optimization and no fma contractions, so that Fortran and C++ results can match bit-for-bit.
# NOTE: we use -O1, not -O0: unoptimized code (especially Kokkos) made the tests 4-15x slower than in
# the regular gnu job, with some tests running out of time. gcc does not reassociate fp operations
# unless asked to (-ffast-math), so -O1 with -ffp-contract=off should not change the results.
SET (ADD_Fortran_FLAGS "-fbacktrace -O1 -ffp-contract=off" CACHE STRING "")
SET (ADD_C_FLAGS       "-O1 -ffp-contract=off" CACHE STRING "")
SET (ADD_CXX_FLAGS     "-O1 -ffp-contract=off" CACHE STRING "")

# The bfb CI jobs only run the cxx_vs_f90 comparisons (and the unit tests), and none of them uses the
# sweqx or tool executables (see .github/workflows/homme-testing.yml), so don't build them.
SET (BUILD_HOMME_SWEQX OFF CACHE BOOL "")
SET (BUILD_HOMME_TOOL  OFF CACHE BOOL "")

include (${CMAKE_CURRENT_LIST_DIR}/ghci-gnu.cmake)
