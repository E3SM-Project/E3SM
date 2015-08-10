cmake_minimum_required (VERSION 3.0.0)
project (PIO C)

#==============================================================================
#  USER-DEFINED OPTIONS (set with "-DOPT=VAL" from command line)
#==============================================================================

option (ENABLE_FORTRAN "Enable the Fortran library builds"  ON)
option (ENABLE_TIMING  "Enable GPTL timing"                 ON)
option (ENABLE_TESTS   "Enable the testing builds"          ON)
option (PREFER_SHARED  "Prefer the use of shared libraries" OFF)

#==============================================================================
#  APPEND TO CMAKE MODULE PATH
#==============================================================================

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================
     
add_subdirectory (src)
add_subdirectory (tests)
