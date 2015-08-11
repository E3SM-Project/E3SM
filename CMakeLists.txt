cmake_minimum_required (VERSION 3.0.0)
project (PIO C)

#==============================================================================
#  USER-DEFINED OPTIONS (set with "-DOPT=VAL" from command line)
#==============================================================================

option (PIO_ENABLE_FORTRAN   "Enable the Fortran library builds"  ON)
option (PIO_ENABLE_TIMING    "Enable GPTL timing"                 ON)
option (PIO_ENABLE_TESTS     "Enable the testing builds"          ON)

option (USE_STATIC_LIBRARIES "Use static libraries for linking"   OFF)

#==============================================================================
#  APPEND TO CMAKE MODULE PATH
#==============================================================================

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================
     
add_subdirectory (src)
add_subdirectory (tests)
