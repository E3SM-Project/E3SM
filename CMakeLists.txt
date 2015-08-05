cmake_minimum_required (VERSION 3.0.0)
project (PIO C)

#==============================================================================
#  USER-DEFINED OPTIONS (set with "-DOPT=VAL" from command line)
#==============================================================================

option (ENABLE_FORTRAN "Enable the Fortran library builds" OFF)
option (ENABLE_TIMING  "Enable GPTL timing in libraries"   OFF)
option (ENABLE_TESTS   "Enable the testing builds"         OFF)

#==============================================================================
#  INCLUDE PATHS
#==============================================================================
     
add_subdirectory (src)
