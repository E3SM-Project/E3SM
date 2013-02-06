PROJECT(HOMME C CXX Fortran)
INCLUDE(FortranCInterface)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# need to get the cmake_fortran_c_interface.h header
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
FortranCInterface_HEADER(cmake_fortran_c_interface.h
                         MACRO_NAMESPACE  "FCI_")
# need to tell the source to use this capability
ADD_DEFINITIONS(-DINCLUDE_CMAKE_FCI)

OPTION(BUILD_HOMME_SWEQX ""   ON)
OPTION(BUILD_HOMME_PREQX ""   ON)
OPTION(BUILD_HOMME_SWDGX ""   ON)
OPTION(BUILD_HOMME_SWIM  ""  OFF)

SET(NUM_POINTS_default   4  CACHE   STRING    "Default number of POINTS (NP) in the element (default value is 4)")
SET(NUM_CELLS_default    4  CACHE   STRING    "Default number of CELLS (NC*NC) in each element (default value is 4)")
SET(NUM_PLEV_default     1  CACHE   STRING    "Default number of POINTS (PLEV) in each element (default value is 1)")
SET(NUM_TRACERS_default  1  CACHE   STRING    "Default number of TRACERS (NTRAC) (default value is 1)")

IF ("${HOMME_ARCH}" STREQUAL "Linux")
  SET(_Linux ON)
  MESSAGE(STATUS "Configuring for Linux")
ELSEIF("${HOMME_ARCH}" STREQUAL "Darwin")
  SET(_Darwin ON)
  MESSAGE(STATUS "Configuring for Darwin")
ENDIF()

IF (NOT DEFINED NUM_POINTS) 
  SET(NUM_POINTS ${NUM_POINTS_default} CACHE STRING "")
ENDIF()
IF (NOT DEFINED NUM_CELLS) 
  SET(NUM_CELLS ${NUM_CELLS_default} CACHE STRING "")
ENDIF()
IF (NOT DEFINED NUM_PLEV) 
  SET(NUM_PLEV ${NUM_PLEV_default} CACHE STRING "")
ENDIF()
IF (NOT DEFINED NUM_TRACERS) 
  SET(NUM_TRACERS ${NUM_TRACERS_default} CACHE STRING "")
ENDIF()
# to allow ifdefs to be included
if(REFSOLN)
  MESSAGE(STATUS "Including option to read or write reference soln data...")
  ADD_DEFINITIONS(-D_REFSOLN)
endif()

MESSAGE(STATUS "Build options:")
MESSAGE(STATUS "-- num points (NP) set to... ${NUM_POINTS}")
MESSAGE(STATUS "-- num cells (NC) set to... ${NUM_CELLS}")
MESSAGE(STATUS "-- num points (PLEV) set to... ${NUM_PLEV}")
MESSAGE(STATUS "-- num tracers (NTRAC) set to... ${NUM_TRACERS}")

# Add the cmake/CMakeLists.txt file
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set up some directories to facilitate finding netcdf and pnetcdf
include(SetHostPaths)

#ADD_SUBDIRECTORY(cmake)

MESSAGE(STATUS "Homme hostname = ${Homme_Hostname}")
MESSAGE(STATUS "Homme OS = ${Homme_OS}")
MESSAGE(STATUS "C Compiler= ${CMAKE_C_COMPILER_ID}")
MESSAGE(STATUS "CXX Compiler= ${CMAKE_CXX_COMPILER_ID}")
MESSAGE(STATUS "Fortran Compiler= ${CMAKE_Fortran_COMPILER_ID}")
MESSAGE(STATUS "CMAKE_SYSTEM= ${CMAKE_SYSTEM}")
MESSAGE(STATUS "CMAKE_SYSTEM_NAME= ${CMAKE_SYSTEM_NAME}")
MESSAGE(STATUS "HOMME_FIND_BLASLAPACK=${HOMME_FIND_BLASLAPACK}")

IF(HOMME_FIND_BLASLAPACK)
  FIND_PACKAGE(BLAS REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
ELSE()
  SET(HOMME_FIND_BLASLAPACK FALSE)
  MESSAGE(STATUS "HOMME_FIND_BLASLAPACK=${HOMME_FIND_BLASLAPACK}")
  # Build the blas and lapack provided with homme
  ADD_SUBDIRECTORY(libs/blas)
  ADD_SUBDIRECTORY(libs/lapack)
ENDIF()

FIND_PACKAGE(Netcdf REQUIRED)
FIND_PACKAGE(Pnetcdf REQUIRED)

ADD_SUBDIRECTORY(utils/pio)
ADD_SUBDIRECTORY(utils/timing)
if(BUILD_HOMME_SWEQX)
  MESSAGE(STATUS "Building sweqx...")
  ADD_SUBDIRECTORY(src/sweqx)
endif()
if(BUILD_HOMME_PREQX)
  MESSAGE(STATUS "Building preqx...")
  ADD_SUBDIRECTORY(src/preqx)
endif()
if(BUILD_HOMME_SWDGX)
  MESSAGE(STATUS "Building swdgx...")
  ADD_SUBDIRECTORY(src/swdgx)
endif()
if(BUILD_HOMME_SWIM)
  MESSAGE(STATUS "Building swim...")
  FIND_PACKAGE(Trilinos REQUIRED)
  MESSAGE(STATUS "Using Trilinos: ${Trilinos_LIBRARY_DIRS}")
  ADD_SUBDIRECTORY(src/swim)
endif()
