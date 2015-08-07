cmake_minimum_required (VERSION 3.0.0)
project (PIO C)

message ("CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")

#==============================================================================
#  USER-DEFINED OPTIONS (set with "-DOPT=VAL" from command line)
#==============================================================================

option (ENABLE_FORTRAN "Enable the Fortran library builds" OFF)
option (ENABLE_TIMING  "Enable GPTL timing"                OFF)
option (ENABLE_TESTS   "Enable the testing builds"         OFF)

#==============================================================================
#  APPEND TO CMAKE MODULE PATH
#==============================================================================

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================
     
add_subdirectory (src)
add_subdirectory (tests)
