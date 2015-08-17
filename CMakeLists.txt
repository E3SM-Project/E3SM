cmake_minimum_required (VERSION 2.8.10)
project (PIO C)

#==============================================================================
#  USER-DEFINED OPTIONS (set with "-DOPT=VAL" from command line)
#==============================================================================

#===== Library Options =====
option (PIO_ENABLE_FORTRAN   "Enable the Fortran library builds"            ON)
option (PIO_ENABLE_TIMING    "Enable GPTL timing"                           ON)
option (PIO_TEST_BIG_ENDIAN  "Enable test to see if machine is big endian"  ON)
option (PIO_ENABLE_MPIIO     "Enable support for MPI-IO auto detect"        ON)

#===== Library Variables
set (PIO_FILESYSTEM_HINTS IGNORE CACHE STRING "Filesystem hints (lustre or gpfs)")

#===== Testing Options =====
option (PIO_ENABLE_TESTS  "Enable the testing builds"                           ON)
option (PIO_HDF5_LOGGING  "Enable hdf5 logging (requires instrumented netcdf4)" OFF)

#==============================================================================
#  APPEND TO CMAKE MODULE PATH
#==============================================================================

#===== Local modules =====
list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#===== External modules =====
if (NOT DEFINED CMAKE_UTILS_DIR)
    message (STATUS "Importing CMake_Fortran_utils")
    execute_process(
        COMMAND git clone https://github.com/CESM-Development/CMake_Fortran_utils
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_QUIET
        ERROR_QUIET)
    find_path (CMAKE_UTILS_DIR
               NAMES mpiexec.cmake
               HINTS ${CMAKE_BINARY_DIR}/CMake_Fortran_utils)
    if (CMAKE_UTILS_DIR)
        message (STATUS "Importing CMake_Fortran_utils - success")
    else ()
        message (FATAL_ERROR "Failed to import CMake_Fortran_utils")
    endif ()
endif ()
set (CMAKE_UTILS_DIR ${CMAKE_UTILS_DIR}
     CACHE STRING "Location of the CMake_Fortran_utils")
list (APPEND CMAKE_MODULE_PATH ${CMAKE_UTILS_DIR})

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================
     
add_subdirectory (src)
if (PIO_ENABLE_TESTS)
    enable_testing()
    include (CTest)
    add_subdirectory (tests)
endif ()
