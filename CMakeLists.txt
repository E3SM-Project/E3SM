PROJECT(HOMME C Fortran)

ENABLE_LANGUAGE(Fortran) 
INCLUDE(FortranCInterface)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

# CTest For Nightly Builds
OPTION(HOMME_USE_CTEST FALSE)
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

# need to tell the source to use this capability
ADD_DEFINITIONS(-DINCLUDE_CMAKE_FCI)

# Set up some directories to facilitate finding netcdf and pnetcdf
INCLUDE(SetHostPaths)
INCLUDE(HommeMacros)

# Default use quad precision
OPTION(HOMME_QUAD_PREC "" TRUE)

IF (${HOMME_QUAD_PREC}) 
  testQuadPrec(HOMME_QUAD_SUPPORTED)
  IF (NOT ${HOMME_QUAD_SUPPORTED})
    SET(HOMME_QUAD_PREC FALSE)
    MESSAGE(STATUS "Quadruple-precision requested but unavailable on this
                    system with this compiler")
  ENDIF()
ENDIF ()

# Debug statement, please leave here
#MESSAGE(STATUS "CMAKE_IMP_LIBS=${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")

# Some options for conditional builds
OPTION(BUILD_HOMME_SWEQX "Shallow water equations FEM" ON)
OPTION(BUILD_HOMME_PREQX "Primitive equations FEM" ON)
OPTION(BUILD_HOMME_SWDGX "Shallow water equations DG" ON)
OPTION(BUILD_HOMME_PRIMDGX "Primitive equations DG" ON)
OPTION(BUILD_HOMME_SWIM  "Shallow water equations implicit" OFF)
OPTION(BUILD_HOMME_PRIM  "Primitive equations implicit" OFF)
# Primitive equations FEM with separate tracer routines
OPTION(BUILD_PREQX_FVM "Primitive equations with FVM tracers" OFF)
OPTION(BUILD_PREQX_SPELT "Primitive euqations with SPELT tracers" OFF)

# Cristoph's "standalone" executables
OPTION(BUILD_HOMME_FVM  "Standalone FVM tracers" OFF)
OPTION(BUILD_HOMME_SPELT  "Standalone SPELT tracers" OFF)

IF (${SWEQX_ONLY})
  SET(BUILD_HOMME_SWEQX "" ON)
  SET(BUILD_HOMME_PREQX "" OFF)
  SET(BUILD_HOMME_PRIM  "" OFF)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_PRIMDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" OFF)
ENDIF ()

IF (${PREQX_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" ON)
  SET(BUILD_HOMME_PRIM  "" OFF)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_PRIMDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" OFF)
ENDIF ()

IF (${SWDGX_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" OFF)
  SET(BUILD_HOMME_PRIM  "" OFF)
  SET(BUILD_HOMME_SWDGX "" ON)
  SET(BUILD_HOMME_PRIMDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" OFF)
ENDIF ()

IF (${PRIMDGX_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" OFF)
  SET(BUILD_HOMME_PRIM  "" OFF)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_PRIMDGX "" ON)
  SET(BUILD_HOMME_SWIM  "" OFF)
ENDIF ()

IF (${SWIM_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" OFF)
  SET(BUILD_HOMME_PRIM  "" OFF)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_PRIMDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" ON)
ENDIF ()

IF (${PRIM_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" OFF)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_PRIMDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" OFF)
  SET(BUILD_HOMME_PRIM  "" ON)
ENDIF ()

IF (${CG_ONLY})
  SET(BUILD_HOMME_SWEQX "" ON)
  SET(BUILD_HOMME_PREQX "" ON)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_PRIMDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" ON)
  SET(BUILD_HOMME_PRIM  "" ON)
ENDIF ()

IF (${DG_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" OFF)
  SET(BUILD_HOMME_SWIM  "" OFF)
  SET(BUILD_HOMME_PRIM  "" OFF)
  SET(BUILD_HOMME_SWDGX "" ON)
  SET(BUILD_HOMME_PRIMDGX "" ON)
ENDIF ()

IF (${3D_CG_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" ON)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" OFF)
  SET(BUILD_HOMME_PRIM  "" ON)
ENDIF ()

OPTION(HOMME_USE_CXX FALSE)
OPTION(HOMME_USE_TRILINOS FALSE)

# swim and prim executables require C++ compiler
IF(${BUILD_HOMME_SWIM} OR ${BUILD_HOMME_PRIM})
  SET (HOMME_USE_CXX TRUE)
  SET (HOMME_USE_TRILINOS TRUE)
ENDIF ()

IF (${HOMME_USE_CXX})
  MESSAGE(STATUS "This configuration of HOMME requires a C++ compiler")
  ENABLE_LANGUAGE(CXX)
ENDIF ()

# Reference solution data
IF(REFSOLN)
  MESSAGE(STATUS "Including option to read or write reference soln data...")
  ADD_DEFINITIONS(-D_REFSOLN)
ENDIF()

MESSAGE(STATUS "Homme hostname = ${Homme_Hostname}")
MESSAGE(STATUS "Homme OS = ${Homme_OS}")
MESSAGE(STATUS "C Compiler = ${CMAKE_C_COMPILER_ID}")
MESSAGE(STATUS "CXX Compiler = ${CMAKE_CXX_COMPILER_ID}")
MESSAGE(STATUS "Fortran Compiler = ${CMAKE_Fortran_COMPILER_ID}")
MESSAGE(STATUS "CMAKE_SYSTEM = ${CMAKE_SYSTEM}")
MESSAGE(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
MESSAGE(STATUS "HOMME_FIND_BLASLAPACK = ${HOMME_FIND_BLASLAPACK}")

OPTION(ENABLE_OPENMP FALSE)
IF (ENABLE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    MESSAGE(STATUS "Found OpenMP Flags")
    # This file is needed for the timing library - this is currently
    # inaccessible from the timing CMake script
    MESSAGE(STATUS "OpenMP_Fortran_FLAGS: ${OpenMP_Fortran_FLAGS}")
    MESSAGE(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
    MESSAGE(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
    MESSAGE(STATUS "OpenMP_EXE_LINKER_FLAGS: ${OpenMP_EXE_LINKER_FLAGS}")
    # The fortran openmp flag should be the same as the C Flag
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  ELSE ()
    MESSAGE(FATAL_ERROR "Unable to find OpenMP")
  ENDIF()
ENDIF ()

#SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${HOMME_ADDITIONAL_FLAGS}")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${HOMME_ADDITIONAL_FLAGS}")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HOMME_ADDITIONAL_FLAGS}")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${HOMME_ADDITIONAL_FLAGS}")

# For PIO Set these options to allow a longer line length (until PIO is fixed)
IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mextend")
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PathScale)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -extend-source")
ENDIF ()

# Intel specific flags
IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O3")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2")
ENDIF ()
 
# Intel specific flags
IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise -ftz")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O3")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2")
 
  # Saving this stuff for later - Do all modern Intel chips support mavx?
  #SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx -DTEMP_INTEL_COMPILER_WORKAROUND_001")
  #SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -DTEMP_INTEL_COMPILER_WORKAROUND_001")
ENDIF ()

# OpenACC specific flag - only supporting PGI compiler
OPTION(ENABLE_OPENACC FALSE)
IF (${ENABLE_OPENACC})
  IF (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
    MESSAGE(FATAL_ERROR "OpenACC only supported through the PGI compiler")
  ELSE ()
    # Need to add -acc to the Fortran FLAGS to see if it will compile 
    # "call acc_init()"
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc")
    TRY_COMPILE(OPENACC_SUCCESS
                ${CMAKE_BINARY_DIR}/compiler_test/
                ${CMAKE_CURRENT_SOURCE_DIR}/cmake/openAccTest.f90
                OUTPUT_VARIABLE COMPILE_OUTPUT)
    IF (${OPENACC_SUCCESS})
      MESSAGE(STATUS "Using OpenACC through PGI compiler")
    ELSE()
      MESSAGE(FATAL_ERROR "Unable to use OpenACC through the PGI compiler")
    ENDIF()
  ENDIF()
ENDIF ()

# CUDA Fortran specific flags - only supporting PGI compiler
OPTION(ENABLE_CUDA_FORTRAN FALSE)
IF (${ENABLE_CUDA_FORTRAN})
  IF (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
    MESSAGE(FATAL_ERROR "CUDA Fortran only supported through the PGI compiler")
  ELSE ()
    # Need to add -Mcuda to the Fortran FLAGS to see if it will compile "use cudafor"
    # Here cc2x is for devices with compute capability >= 2.0 
    # To do: determine a way to generalize this
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ta=nvidia -Mcuda=4.2,cc2x,ptxinfo,keepgpu")
    TRY_COMPILE(CUDAFOR
                ${CMAKE_BINARY_DIR}/compiler_test/
                ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cudaFortranTest.F90
                OUTPUT_VARIABLE COMPILE_OUTPUT)
    IF (${CUDAFOR})
      SET(PREQX_USE_CUDA_FORTRAN TRUE)
      MESSAGE(STATUS "Using CUDA Fortran through PGI compiler")
    ELSE()
      SET(PREQX_USE_CUDA_FORTRAN FALSE)
      MESSAGE(FATAL_ERROR "Unable to use CUDA Fortran through the PGI compiler")
    ENDIF()
  ENDIF()
ENDIF ()

# Allow the option to add compiler flags to those provided
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ADD_Fortran_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ADD_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADD_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADD_LINKER_FLAGS}")

# Allow the option to override compiler flags
IF (FORCE_Fortran_FLAGS)
  SET(CMAKE_Fortran_FLAGS ${FORCE_Fortran_FLAGS})
ENDIF ()
IF (FORCE_C_FLAGS)
  SET(CMAKE_C_FLAGS ${FORCE_C_FLAGS})
ENDIF ()
IF (FORCE_CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${FORCE_CXX_FLAGS})
ENDIF ()
IF (FORCE_LINKER_FLAGS)
  SET(CMAKE_EXE_LINKER_FLAGS ${FORCE_LINKER_FLAGS})
ENDIF ()

# Print the compiler flags to be used
MESSAGE(STATUS "Fortran Flags = ${CMAKE_Fortran_FLAGS}")
MESSAGE(STATUS "C Flags = ${CMAKE_C_FLAGS}")
MESSAGE(STATUS "CXX Flags = ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Linker Flags = ${CMAKE_EXE_LINKER_FLAGS}")

OPTION(HOMME_FIND_BLASLAPACK FALSE)

IF(${HOMME_FIND_BLASLAPACK})
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
OPTION(WITH_PNETCDF TRUE)
IF (${WITH_PNETCDF})
  FIND_PACKAGE(Pnetcdf REQUIRED)
ENDIF ()

IF (${HOMME_USE_TRILINOS})
  MESSAGE(STATUS "This configuration of requires Trilinos")
  FIND_PACKAGE(Trilinos REQUIRED)
ENDIF ()

##############################################################################
# Compiler test to see if the mpi compiler supports "use mpi"
# Note 1) Although "use mpi" is part of the standard, many installs of MPI are
#   incorrect and this doesn't work.
# Note 2) Steve G. has pushed this exact fix into PIO, which is the only place
#   that needs it. Once that code is in place this code will be deleted.
TRY_COMPILE(MPIMOD_SUCCESS 
            ${CMAKE_BINARY_DIR}/compiler_test/ ${CMAKE_CURRENT_SOURCE_DIR}/cmake/mpiModTest.f90
            OUTPUT_VARIABLE COMPILE_OUTPUT)
IF (${MPIMOD_SUCCESS})
  SET(PIO_USE_MPIMOD TRUE)
  MESSAGE(STATUS "MPI Fortran \"use\" mpi module available.")
ELSE()
  SET(PIO_USE_MPIMOD FALSE)
  MESSAGE(STATUS "MPI Fortran \"use\" mpi module not available.")
ENDIF()
##############################################################################
IF (NOT PIO_USE_MPIMOD)
  ADD_DEFINITIONS(-DNO_MPIMOD)
ENDIF()
##############################################################################

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


