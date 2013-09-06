
PROJECT(HOMME C Fortran)

ENABLE_LANGUAGE(Fortran) 
INCLUDE(FortranCInterface)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

# CTest For Nightly/Experimental Builds
OPTION(HOMME_USE_CTEST "Use CTest Dashboard" FALSE)
IF (${HOMME_USE_CTEST})
  INCLUDE(CTest)
ENDIF ()

# For testing
ENABLE_TESTING()

# Add the cmake/CMakeLists.txt file
SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# need to get the cmake_fortran_c_interface.h header
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
FortranCInterface_HEADER(cmake_fortran_c_interface.h
                         MACRO_NAMESPACE  "FCI_")

# Prefer shared libraries (Yellowstone etc.) 
OPTION(PREFER_SHARED "Prefer linking with shared libraries." FALSE)

# need to tell the source to use this capability
ADD_DEFINITIONS(-DINCLUDE_CMAKE_FCI)

# Set up some directories to facilitate finding netcdf and pnetcdf
INCLUDE(SetHostPaths)
INCLUDE(HommeMacros)

# Debug statement, please leave here
#MESSAGE(STATUS "CMAKE_IMP_LIBS=${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")

##############################################################################
# Compiler test to see if the mpi compiler supports "use mpi" or "include
# mpi.h"
# Note 1) Although "use mpi" is part of the standard, many installs of MPI are
#   incorrect and this doesn't work.
# Note 2) Steve G. has pushed this exact fix into PIO, which is the only place
#   that needs it. Once that code is in place this code will be deleted.
##############################################################################
TRY_COMPILE(MPIMOD_SUCCESS 
            ${CMAKE_BINARY_DIR}/tests/compilerTests/
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilerTests/mpiModTest.f90
            OUTPUT_VARIABLE COMPILE_OUTPUT)
IF (${MPIMOD_SUCCESS})
  SET(PIO_USE_MPIMOD TRUE)
  MESSAGE(STATUS "MPI Fortran \"use\" mpi module available.")
ELSE()
  SET(PIO_USE_MPIMOD FALSE)
  MESSAGE(STATUS "MPI Fortran \"use\" mpi module not available.")
  TRY_COMPILE(MPIINC_SUCCESS 
            ${CMAKE_BINARY_DIR}/tests/compilerTests/
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilerTests/mpiIncludeTest.F90
            OUTPUT_VARIABLE COMPILE_OUTPUT)
  IF(NOT ${MPIINC_SUCCESS})
    MESSAGE(FATAL_ERROR "MPI fortran \"include mpi.h\" not available.\n"
            "Fortran Compiler not capable of compiling MPI programs.\n"
            "Be sure to use a parallel (MPI) compiler.")
  ENDIF ()
ENDIF()
##############################################################################
IF (NOT PIO_USE_MPIMOD)
  ADD_DEFINITIONS(-DNO_MPIMOD)
ENDIF()
##############################################################################

# Default use quad precision
OPTION(HOMME_QUAD_PREC "Default to try to use quadruple precision" TRUE)

IF (${HOMME_QUAD_PREC}) 
  testQuadPrec(HOMME_QUAD_SUPPORTED)
  IF (NOT ${HOMME_QUAD_SUPPORTED})
    SET(HOMME_QUAD_PREC FALSE)
    MESSAGE(STATUS "Quadruple-precision requested but unavailable on this
                    system with this compiler")
  ENDIF()
ENDIF ()

# Some options for conditional builds
OPTION(BUILD_HOMME_SWEQX "Shallow water equations FEM" ON)
OPTION(BUILD_HOMME_PREQX "Primitive equations FEM" ON)
OPTION(BUILD_HOMME_SWDGX "Shallow water equations DG" ON)
OPTION(BUILD_HOMME_PRIMDGX "Primitive equations DG" ON)
OPTION(BUILD_HOMME_SWIM  "Shallow water equations implicit" OFF)
OPTION(BUILD_HOMME_PRIM  "Primitive equations implicit" OFF)

# Primitive equations FEM with separate tracer routines
OPTION(BUILD_PREQX_SPELT "Primitive euqations with SPELT tracers" OFF)

# Cristoph's "standalone" executables
OPTION(BUILD_HOMME_FVM  "Standalone FVM tracers" OFF)
OPTION(BUILD_HOMME_SPELT  "Standalone SPELT tracers" OFF)

# swim and prim executables require C++ compiler
IF(${BUILD_HOMME_SWIM} OR ${BUILD_HOMME_PRIM})
  SET (HOMME_USE_CXX TRUE)
  SET (HOMME_USE_TRILINOS TRUE)
ELSE ()
  SET(HOMME_USE_CXX FALSE)
  SET(HOMME_USE_TRILINOS FALSE)
ENDIF ()

IF (${HOMME_USE_CXX})
  MESSAGE(STATUS "This configuration of HOMME requires a C++ compiler")
  ENABLE_LANGUAGE(CXX)
ENDIF ()

OPTION(REFSOLN "Whether to read/write reference soln data" FALSE)
# Reference solution data
IF(REFSOLN)
  MESSAGE(STATUS "Including option to read or write reference soln data...")
  ADD_DEFINITIONS(-D_REFSOLN)
ENDIF()

##############################################################################
# Compiler options set in cmake/SetCompiler.cmake
##############################################################################
INCLUDE(SetCompilerFlags)

# Print the compiler flags to be used
MESSAGE(STATUS "Fortran Flags = ${CMAKE_Fortran_FLAGS}")
MESSAGE(STATUS "C Flags = ${CMAKE_C_FLAGS}")
MESSAGE(STATUS "CXX Flags = ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Linker Flags = ${CMAKE_EXE_LINKER_FLAGS}")
##############################################################################

OPTION(HOMME_FIND_BLASLAPACK "Whether to use system blas/lapack" FALSE)
MESSAGE(STATUS "HOMME_FIND_BLASLAPACK=${HOMME_FIND_BLASLAPACK}")

IF(${HOMME_FIND_BLASLAPACK})
  FIND_PACKAGE(BLAS REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
ELSE()
  SET(HOMME_FIND_BLASLAPACK FALSE)
  # Build the blas and lapack provided with homme
  ADD_SUBDIRECTORY(libs/blas)
  ADD_SUBDIRECTORY(libs/lapack)
ENDIF()

# Netcdf is required
FIND_PACKAGE(Netcdf REQUIRED)

# PNetcdf is optional but used by default
OPTION(WITH_PNETCDF "Whether to build with PNetcdf" TRUE)
IF (${WITH_PNETCDF})
  FIND_PACKAGE(Pnetcdf REQUIRED)
ELSE ()
  MESSAGE(WARNING "Warning: Not building with PNetcdf - cannot run all regression tests.")
ENDIF ()

IF (${HOMME_USE_TRILINOS})
  MESSAGE(STATUS "This configuration of requires Trilinos")
  FIND_PACKAGE(Trilinos REQUIRED)
ENDIF ()

# Zoltan partitioning library
IF (DEFINED ZOLTAN_DIR)
  MESSAGE(STATUS "Building with Zoltan")
  FIND_PACKAGE(Zoltan REQUIRED)
  SET(HAVE_ZOLTAN TRUE)
ENDIF ()

# This is required in timing to set HAVE_MPI
ADD_DEFINITIONS(-DSPMD)

SET(PIO_USE_MPIIO OFF CACHE BOOL "Turn this off due to bug in pio")

ADD_SUBDIRECTORY(utils/pio)

ADD_SUBDIRECTORY(utils/timing)

IF(${BUILD_HOMME_SWEQX})
  ADD_SUBDIRECTORY(src/sweqx)
ENDIF()

IF(${BUILD_HOMME_PREQX})
  ADD_SUBDIRECTORY(src/preqx)
ENDIF()

IF(${BUILD_HOMME_SWDGX})
  ADD_SUBDIRECTORY(src/swdgx)
ENDIF()

IF(${BUILD_HOMME_PRIMDGX})
  ADD_SUBDIRECTORY(src/primdgx)
ENDIF()

IF(${BUILD_HOMME_FVM})
  ADD_SUBDIRECTORY(src/fvm)
ENDIF()

IF(${BUILD_HOMME_SPELT})
  ADD_SUBDIRECTORY(src/spelt)
ENDIF()

IF(${BUILD_HOMME_SWIM})
  ADD_SUBDIRECTORY(src/swim)
ENDIF()

IF(${BUILD_HOMME_PRIM})
  ADD_SUBDIRECTORY(src/prim)
ENDIF()

ADD_SUBDIRECTORY(test_execs)


