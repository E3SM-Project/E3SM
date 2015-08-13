cmake_minimum_required (VERSION 3.0.0)
project (PIO C)

#==============================================================================
#  USER-DEFINED OPTIONS (set with "-DOPT=VAL" from command line)
#==============================================================================

#===== Library Options =====
option (PIO_ENABLE_FORTRAN    "Enable the Fortran library builds"    ON)
option (PIO_ENABLE_TIMING     "Enable GPTL timing"                   ON)

#===== Testing Options =====
option (PIO_ENABLE_TESTS  "Enable the testing builds"                           ON)
option (PIO_HDF5_LOGGING  "Enable hdf5 logging (requires instrumented netcdf4)" OFF)

#==============================================================================
#  APPEND TO CMAKE MODULE PATH
#==============================================================================

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================
     
add_subdirectory (src)
if (PIO_ENABLE_TESTS)
    add_subdirectory (tests)
endif ()
