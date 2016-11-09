
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

OPTION(IS_BIG_ENDIAN "Whether the machine is big endian" FALSE)
#TRY_RUN(IS_BIG_ENDIAN BIG_ENDIAN_COMPILED
#  ${CMAKE_BINARY_DIR}/tests/compilerTests/
#  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilerTests/isBigEndian.c)

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
OPTION(BUILD_HOMME_PREQX_ACC "Primitive equations FEM with OpenACC" ON)
OPTION(BUILD_HOMME_SWIM  "Shallow water equations implicit" OFF)
OPTION(BUILD_HOMME_PRIM  "Primitive equations implicit" OFF)

# Cristoph's "standalone" executables
OPTION(BUILD_HOMME_FVM  "Standalone FVM tracers" OFF)

# swim and prim executables require C++ compiler
IF(${BUILD_HOMME_SWIM} OR ${BUILD_HOMME_PRIM})
  SET (HOMME_USE_CXX TRUE)
  SET (HOMME_USE_TRILINOS TRUE)
ELSE ()
  SET(HOMME_USE_CXX FALSE)
ENDIF ()

IF (${HOMME_USE_TRILINOS}) 
  SET (HOMME_USE_CXX TRUE)
ENDIF()


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


OPTION(ENABLE_PERSISTENT_MPI "Enable persistent MPI" OFF)
IF (${ENABLE_PERSISTENT_MPI})
  SET(MPI_PERSISTENT TRUE)
ENDIF ()

OPTION(HOMME_USE_MKL "Whether to use Intel's MKL instead of blas/lapack" FALSE)
MESSAGE(STATUS "HOMME_USE_MKL=${HOMME_USE_MKL}")
IF(${HOMME_USE_MKL})
  MESSAGE(STATUS "For Edison/Cori at NERSC, MKL included with Intel compiler")
ELSE()
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
ENDIF()

IF (${HOMME_USE_TRILINOS})
  MESSAGE(STATUS "This configuration of requires Trilinos")
  FIND_PACKAGE(Trilinos REQUIRED)
  SET(HAVE_TRILINOS TRUE)
  IF("${Trilinos_PACKAGE_LIST}"  MATCHES ".*Zoltan2*.")
    MESSAGE(STATUS "Trilinos is compiled with Zoltan2 Trilinos_PACKAGE_LIST:${Trilinos_PACKAGE_LIST}")
    MESSAGE(STATUS "ENABLING ZOLTAN2")
    SET(TRILINOS_HAVE_ZOLTAN2 TRUE)
  ENDIF ()
  add_definitions(${Trilinos_CXX_COMPILER_FLAGS})
  include_directories(${Trilinos_TPL_INCLUDE_DIRS})
  include_directories(${Trilinos_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})
  link_directories(${Trilinos_TPL_LIBRARY_DIRS})
ENDIF ()

# Print the compiler flags to be used
MESSAGE(STATUS "\nFortran Flags = ${CMAKE_Fortran_FLAGS}")
MESSAGE(STATUS "C Flags = ${CMAKE_C_FLAGS}")
MESSAGE(STATUS "CXX Flags = ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Linker Flags = ${CMAKE_EXE_LINKER_FLAGS} \n")
##############################################################################


# Zoltan partitioning library
IF (DEFINED ZOLTAN_DIR)
  MESSAGE(STATUS "Building with Zoltan")
  FIND_PACKAGE(Zoltan REQUIRED)
  SET(HAVE_ZOLTAN TRUE)
ENDIF ()

IF (DEFINED EXTRAE_DIR)
  MESSAGE(STATUS "Building with Extrae")
  FIND_PACKAGE(Extrae REQUIRED)
  SET(HAVE_EXTRAE TRUE)
ENDIF ()

###########################################
# Definitions for timing library
###########################################
# This is required in timing to set HAVE_MPI
ADD_DEFINITIONS(-DSPMD)
OPTION(ENABLE_NANOTIMERS "Use nano timers in timing library" FALSE)
OPTION(USE_BIT64 "Set BIT64 (for 64 bit arch) in timing library when ENABLE_NANOTIME" FALSE)
IF (${ENABLE_NANOTIMERS})
  ADD_DEFINITIONS(-DHAVE_NANOTIME)
  IF (${USE_BIT64})
    ADD_DEFINITIONS(-DBIT64)
  ENDIF ()
ENDIF ()

# Use ccsm_utils FindNETCDF path
SET (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../cime/externals/CMake" ${CMAKE_MODULE_PATH})
FIND_PACKAGE(NETCDF)
# For pio - allow the user to set HDF5_DIR to find libhdf5 etc.
IF (HDF5_DIR)
  SET(ENV{HDF5_ROOT} ${HDF5_DIR})
ENDIF ()
ADD_SUBDIRECTORY(utils/pio)

# CPRNC locations
SET (CPRNC_INSTALL_DIR ${CMAKE_BINARY_DIR}/utils/cprnc)
SET (CPRNC_BINARY ${CMAKE_BINARY_DIR}/utils/cprnc/cprnc)
ADD_SUBDIRECTORY(utils/cprnc)

ADD_SUBDIRECTORY(utils/timing)

SET(EXEC_LIST ""  CACHE INTERNAL "List of configured executables")
IF(${BUILD_HOMME_SWEQX})
  ADD_SUBDIRECTORY(src/sweqx)
ENDIF()

IF(${BUILD_HOMME_PREQX})
  ADD_SUBDIRECTORY(src/preqx)
ENDIF()

IF(${BUILD_HOMME_PREQX_ACC})
  ADD_SUBDIRECTORY(src/preqx_acc)
ENDIF()

IF(${BUILD_HOMME_FVM})
  ADD_SUBDIRECTORY(src/fvm)
ENDIF()

IF(${BUILD_HOMME_SWIM})
  ADD_SUBDIRECTORY(src/swim)
ENDIF()

IF(${BUILD_HOMME_PRIM})
  ADD_SUBDIRECTORY(src/prim)
ENDIF()

ADD_SUBDIRECTORY(test_execs)

# Create a "make clean-all" target to delete CMake cache and the cprnc build
ADD_CUSTOM_TARGET(clean-all
   COMMAND ${CMAKE_BUILD_TOOL} clean
   COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/clean-cmake-files.cmake"
)

ADD_SUBDIRECTORY(dcmip_tests)
