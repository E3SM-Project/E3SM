cmake_minimum_required (VERSION 2.8.12)
project (PIO C Fortran)
#cmake_policy(VERSION 3.5.2)

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
option (PIO_ENABLE_LOGGING   "Enable debug logging (large output possible)" OFF)
option (PIO_ENABLE_DOC       "Enable building PIO documentation"            ON)
option (PIO_ENABLE_COVERAGE  "Enable code coverage"                         OFF)
option (PIO_ENABLE_EXAMPLES  "Enable PIO examples"                          ON)
option (PIO_INTERNAL_DOC     "Enable PIO developer documentation"           OFF)
option (PIO_TEST_BIG_ENDIAN  "Enable test to see if machine is big endian"  ON)
option (PIO_USE_MPIIO        "Enable support for MPI-IO auto detect"        ON)
option (PIO_USE_MPISERIAL    "Enable mpi-serial support (instead of MPI)"   OFF)
option (PIO_USE_MALLOC       "Use native malloc (instead of bget package)"  OFF)
option (WITH_PNETCDF         "Require the use of PnetCDF"                   ON)

# Set a variable that appears in the config.h.in file.
if(PIO_USE_MALLOC)
  set(USE_MALLOC 1)
else()
  set(USE_MALLOC 0)
endif()

# Set a variable that appears in the config.h.in file.
if(PIO_ENABLE_LOGGING)
  set(ENABLE_LOGGING 1)
else()
  set(ENABLE_LOGGING 0)
endif()

if(PIO_USE_MPISERIAL)
  set(USE_MPI_SERIAL 1)
else()
  set(USE_MPI_SERIAL 0)
endif()

#===== Library Variables =====
set (PIO_FILESYSTEM_HINTS IGNORE CACHE STRING "Filesystem hints (lustre or gpfs)")

#===== Testing Options =====
option (PIO_ENABLE_TESTS  "Enable the testing builds"                           ON)
option (PIO_VALGRIND_CHECK  "Enable memory leak check using valgrind"           OFF)

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
#  SET CODE COVERAGE COMPILER FLAGS
#==============================================================================

# Only support GNU compilers at this time
if (PIO_ENABLE_COVERAGE)
  if (CMAKE_C_COMPILER_NAME STREQUAL "GNU")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
  else ()
    message (WARNING "The C compiler is non-GNU: coverage of C code could NOT be enabled")
  endif ()
  if (CMAKE_Fortran_COMPILER_NAME STREQUAL "GNU")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage")
  else ()
    message (WARNING "The Fortran compiler is non-GNU: coverage of Fortran code could NOT be enabled")
  endif ()
endif ()

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================

# Libraries
add_subdirectory (src)

#==============================================================================
#  TESTING TARGET
#==============================================================================

# Custom "piotests" target (builds the test executables)
add_custom_target (tests)

# Custom "check" target that depends upon "tests"
add_custom_target (check COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies (check tests)

# Tests
if (PIO_ENABLE_TESTS)
  enable_testing()
  include (CTest)
  add_subdirectory (tests)
  if (PIO_ENABLE_EXAMPLES)
    add_subdirectory (examples)
  endif ()
endif ()

# Documentation
if (PIO_ENABLE_DOC)
  add_subdirectory (doc)
endif ()
