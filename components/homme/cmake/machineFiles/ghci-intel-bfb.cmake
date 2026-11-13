# CMake initial cache file for the ghci intel image, for BFB (cxx vs f90) testing
#
# NOTE: these must be set BEFORE including the non-bfb file (the first value set in the cache wins)
SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

# Strict fp model and low optimization, so that Fortran and C++ results can match bit-for-bit
SET (ADD_Fortran_FLAGS "-traceback -fp-model strict -O1" CACHE STRING "")
SET (ADD_C_FLAGS       "-fp-model strict -O1" CACHE STRING "")
# Also turn off annoying compiler warning (same as in the non-bfb file)
SET (ADD_CXX_FLAGS     "-fp-model strict -O1 -Wno-overriding-option" CACHE STRING "")

# The bfb CI jobs only run the cxx_vs_f90 comparisons (and the unit tests), and none of them uses the
# sweqx or tool executables (see .github/workflows/homme-testing.yml), so don't build them.
SET (BUILD_HOMME_SWEQX OFF CACHE BOOL "")
SET (BUILD_HOMME_TOOL  OFF CACHE BOOL "")

include (${CMAKE_CURRENT_LIST_DIR}/ghci-intel.cmake)
