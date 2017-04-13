cmake_minimum_required (VERSION 2.8.12)
project (PIO C Fortran)

# The project version number.
set(VERSION_MAJOR   2   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(VERSION_PATCH  28   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

#==============================================================================
#  USER-DEFINED OPTIONS (set with "-DOPT=VAL" from command line)
#==============================================================================

#===== Library Options =====
option (PIO_ENABLE_FORTRAN   "Enable the Fortran library builds"            ON)
option (PIO_ENABLE_TIMING    "Enable the use of the GPTL timing library"    ON)
option (PIO_TEST_BIG_ENDIAN  "Enable test to see if machine is big endian"  ON)
option (PIO_USE_MPIIO        "Enable support for MPI-IO auto detect"        ON)
option (PIO_USE_MPISERIAL    "Enable mpi-serial support (instead of MPI)"   OFF)
option (PIO_USE_MALLOC       "Use native malloc (instead of bget package)"  OFF)
option (WITH_PNETCDF         "Require the use of PnetCDF"                   ON)

if(PIO_USE_MALLOC)
  set(USE_MALLOC 1)
else()
  set(USE_MALLOC 0)
endif()

#===== Library Variables =====
set (PIO_FILESYSTEM_HINTS IGNORE CACHE STRING "Filesystem hints (lustre or gpfs)")

#===== Testing Options =====
option (PIO_ENABLE_TESTS  "Enable the testing builds"                           ON)
option (PIO_HDF5_LOGGING  "Enable hdf5 logging (requires instrumented netcdf4)" OFF)

#==============================================================================
#  BACKWARDS COMPATIBILITY
#==============================================================================

# Old NETCDF_DIR variable --> NetCDF_PATH
if (DEFINED NETCDF_DIR)
    set (NetCDF_PATH ${NETCDF_DIR}
         CACHE STRING "Location of the NetCDF library installation")
endif ()

# Old PNETCDF_DIR variable --> PnetCDF_PATH
if (DEFINED PNETCDF_DIR)
    set (PnetCDF_PATH ${PNETCDF_DIR}
         CACHE STRING "Location of the PnetCDF library installation")
endif ()

#==============================================================================
#  PREPEND TO CMAKE MODULE PATH
#==============================================================================

#===== Local modules =====
list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#===== External modules =====
if (NOT DEFINED USER_CMAKE_MODULE_PATH)
    message (STATUS "Importing CMake_Fortran_utils")
    execute_process(
        COMMAND git clone https://github.com/CESM-Development/CMake_Fortran_utils
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_QUIET
        ERROR_QUIET)
    find_path (USER_CMAKE_MODULE_PATH
               NAMES mpiexec.cmake
               HINTS ${CMAKE_BINARY_DIR}/CMake_Fortran_utils)
    if (USER_CMAKE_MODULE_PATH)
        message (STATUS "Importing CMake_Fortran_utils - success")
    else ()
        message (FATAL_ERROR "Failed to import CMake_Fortran_utils")
    endif ()
endif ()
set (USER_CMAKE_MODULE_PATH ${USER_CMAKE_MODULE_PATH}
     CACHE STRING "Location of the CMake_Fortran_utils")
list (APPEND CMAKE_MODULE_PATH ${USER_CMAKE_MODULE_PATH})

#==============================================================================
#  HELPFUL GLOBAL VARIABLES
#==============================================================================

# System Name
string (TOUPPER "${CMAKE_SYSTEM_NAME}" CMAKE_SYSTEM_NAME_CAPS)
set (CMAKE_SYSTEM_DIRECTIVE "${CMAKE_SYSTEM_NAME_CAPS}"
     CACHE STRING "System name preprocessor directive")

# C Compiler Name
string (TOUPPER "${CMAKE_C_COMPILER_ID}" CMAKE_C_COMPILER_NAME)
if (CMAKE_C_COMPILER_NAME STREQUAL "XL")
    set (CMAKE_C_COMPILER_NAME "IBM")
endif ()
set (CMAKE_C_COMPILER_DIRECTIVE "CPR${CMAKE_C_COMPILER_NAME}"
     CACHE STRING "C compiler name preprocessor directive")

# Fortran Compiler Name
string (TOUPPER "${CMAKE_Fortran_COMPILER_ID}" CMAKE_Fortran_COMPILER_NAME)
if (CMAKE_Fortran_COMPILER_NAME STREQUAL "XL")
    set (CMAKE_Fortran_COMPILER_NAME "IBM")
endif ()
set (CMAKE_Fortran_COMPILER_DIRECTIVE "CPR${CMAKE_Fortran_COMPILER_NAME}"
     CACHE STRING "Fortran compiler name preprocessor directive")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/clib/config.h.in"
  "${PROJECT_BINARY_DIR}/src/clib/config.h"
  )

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================

# Libraries
add_subdirectory (src)

# Tests
if (PIO_ENABLE_TESTS)
    enable_testing()
    include (CTest)
    add_subdirectory (tests)
    add_subdirectory (examples)
endif ()

# Documentation
add_subdirectory (doc)
