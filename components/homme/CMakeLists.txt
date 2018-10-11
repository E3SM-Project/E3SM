# Best practice to put cmake required version at the top
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

PROJECT(HOMME C Fortran)

ENABLE_LANGUAGE(Fortran) 

# Print the sha of the last commit,
# useful to double check which version was build in the build folder.
EXECUTE_PROCESS (COMMAND git rev-parse HEAD
                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                   OUTPUT_VARIABLE LAST_GIT_COMMIT_SHA
                   OUTPUT_STRIP_TRAILING_WHITESPACE)
SET (LAST_GIT_COMMIT_SHA ${LAST_GIT_COMMIT_SHA} CACHE STRING "The sha of the last git commit.")
MESSAGE (STATUS "The sha of the last commit is ${LAST_GIT_COMMIT_SHA}")

INCLUDE(FortranCInterface)

# CTest For Nightly/Experimental Builds
OPTION(HOMME_USE_CTEST "Use CTest Dashboard" FALSE)
IF (${HOMME_USE_CTEST})
  INCLUDE(CTest)
ENDIF ()

# For testing
ENABLE_TESTING()
SET (HOMME_TESTING_PROFILE "nightly" CACHE STRING "Determine how long/pervasive the testing is. Currently available options: 'dev', 'short', 'nightly'.")

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

# In Kokkos builds, tell Config.hpp that CMake is being used to build.
add_definitions(-DHOMMEXX_CONFIG_IS_CMAKE)

SET (PERFORMANCE_PROFILE "" CACHE STRING "Whether to build and link with various profiler libraries")

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
OPTION(BUILD_HOMME_THETA "SB81 NH model with theta and Exner pressure" ON)
OPTION(BUILD_HOMME_PREQX_ACC "Primitive equations FEM with OpenACC" ON)
OPTION(BUILD_HOMME_PREQX_KOKKOS "Primitive equations FEM with Kokkos" OFF)
OPTION(BUILD_HOMME_PESE  "Primitive equations with vertical spectral elements" OFF)
OPTION(BUILD_HOMME_SWIM  "Shallow water equations implicit" OFF)
OPTION(BUILD_HOMME_PRIM  "Primitive equations implicit" OFF)

#by default we don't need cxx
SET(HOMME_USE_CXX FALSE)

# swim and prim executables require C++ compiler
IF(${BUILD_HOMME_SWIM} OR ${BUILD_HOMME_PRIM})
  SET (HOMME_USE_CXX TRUE)
  SET (HOMME_USE_TRILINOS TRUE)
ENDIF ()
#if trilinos is needed, we need cxx
IF (${HOMME_USE_TRILINOS})
  SET (HOMME_USE_CXX TRUE)
ENDIF()
#if kokkos is needed, then so is cxx
SET (HOMME_USE_KOKKOS FALSE)
IF (${BUILD_HOMME_PREQX_KOKKOS})
  SET (HOMME_USE_KOKKOS TRUE)
  SET (HOMME_USE_CXX TRUE)
  IF (NOT ${BUILD_HOMME_PREQX})
    # If we build theta kokkos, we also build theta, so we can compare
    MESSAGE (STATUS "Setting manually disabled BUILD_HOMME_PREQX to ON, since BUILD_HOMME_PREQX_KOKKOS is ON.")
  ENDIF ()

  # We need to compare against F90 implementation, so turn this on
  SET (BUILD_HOMME_PREQX TRUE)
ENDIF ()

IF (${HOMME_USE_CXX})
  MESSAGE(STATUS "This configuration of HOMME requires a C++ compiler")
  ENABLE_LANGUAGE(CXX)
ENDIF ()

IF (${HOMME_USE_KOKKOS})
  ############################################
  # Selection of Kokkos execution space
  ############################################

  SET (HOMMEXX_EXEC_SPACE "Default" CACHE STRING "Select the kokkos exec space")

  STRING (TOUPPER ${HOMMEXX_EXEC_SPACE} HOMMEXX_EXEC_SPACE_UPPER)
  IF (${HOMMEXX_EXEC_SPACE_UPPER} STREQUAL "CUDA")
    SET (HOMMEXX_CUDA_SPACE ON)
  ELSEIF (${HOMMEXX_EXEC_SPACE_UPPER} STREQUAL "OPENMP")
    SET (HOMMEXX_OPENMP_SPACE ON)
  ELSEIF (${HOMMEXX_EXEC_SPACE_UPPER} STREQUAL "THREADS")
    SET (HOMMEXX_THREADS_SPACE ON)
  ELSEIF (${HOMMEXX_EXEC_SPACE_UPPER} STREQUAL "SERIAL")
    SET (HOMMEXX_SERIAL_SPACE ON)
  ELSEIF (${HOMMEXX_EXEC_SPACE_UPPER} STREQUAL "DEFAULT")
    SET (HOMMEXX_DEFAULT_SPACE ON)
  ELSE()
    MESSAGE (ABORT "Invalid choice for 'HOMMEXX_EXEC_SPACE'. Valid options (case insensitive) are 'Cuda', 'OpenMP', 'Threads', 'Serial', 'Default'")
  ENDIF()

  # Execution space parameters. 8 is a nice size for V100.
  SET (HOMMEXX_CUDA_MIN_WARP_PER_TEAM 8 CACHE INT "Minimum number of warps to get 100% occoupancy on GPU")
  SET (HOMMEXX_CUDA_MAX_WARP_PER_TEAM 16 CACHE INT "Maximum number of warps to get 100% occoupancy on GPU")

  # An option to allow to use GPU pointers for MPI calls. The value of this option is irrelevant for CPU/KNL builds.
  OPTION (HOMMEXX_MPI_ON_DEVICE "Whether we want to use device pointers for MPI calls (relevant only for GPU builds)" ON)
ENDIF()

##############################################################################
# Compiler options set in cmake/SetCompiler.cmake
##############################################################################
INCLUDE(SetCompilerFlags)


OPTION(HOMME_USE_MKL "Whether to use Intel's MKL instead of blas/lapack" FALSE)
IF(${HOMME_USE_MKL})
  MESSAGE(STATUS "HOMME_USE_MKL=${HOMME_USE_MKL}, adding -mkl")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mkl")
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

MESSAGE(STATUS "Fortran Flags = ${CMAKE_Fortran_FLAGS}")
MESSAGE(STATUS "C Flags = ${CMAKE_C_FLAGS}")
MESSAGE(STATUS "CXX Flags = ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Linker Flags = ${CMAKE_EXE_LINKER_FLAGS}")

IF (${HOMME_USE_KOKKOS})
  IF (${CUDA_FOUND})
    # Note: If you ask for GPU to be bfb with CPU, you should also ask for a strict fp model on cpu.
    OPTION (HOMMEXX_GPU_BFB_WITH_CPU "Whether we want cpu and gpu solution to be bfb" OFF)
  ENDIF ()

  SET (HOMMEXX_VECTOR_SIZE 8 CACHE STRING
    "If AVX or Cuda don't take priority, use this software vector size.")

  IF (CMAKE_BUILD_TYPE_UPPER MATCHES "DEBUG" OR CMAKE_BUILD_TYPE_UPPER MATCHES "RELWITHDEBINFO")
    SET (HOMMEXX_DEBUG ON)
  ENDIF()

  CONFIGURE_FILE (${CMAKE_CURRENT_SOURCE_DIR}/src/share/cxx/Hommexx_config.h.in ${CMAKE_BINARY_DIR}/src/share/cxx/Hommexx_config.h)

ENDIF ()


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
ADD_SUBDIRECTORY(utils/cime/src/share/timing)


# CMAKE_CURRENT_SOURCE_DIR == homme
SET (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
SET (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/utils/cime/src/CMake" ${CMAKE_MODULE_PATH})
SET (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/utils/cime/src/externals/pio2/cmake" ${CMAKE_MODULE_PATH})



#
# PIO
#
#set search paths for PIO's findNetCDF
SET(NetCDF_PATH ${NETCDF_DIR})
SET(PnetCDF_PATH ${PNETCDF_DIR})
message("** Configuring PIO with:")
message("-- NetCDF_PATH = ${NetCDF_PATH}")
message("-- PnetCDF_PATH = ${PnetCDF_PATH}")
# pio needs cime/externals/genf90/genf90.pl
SET(GENF90_PATH ${CMAKE_CURRENT_SOURCE_DIR}/utils/cime/src/externals/genf90)
ADD_SUBDIRECTORY(utils/cime/src/externals/pio1)



# CPRNC locations
IF (CPRNC_DIR)
  # location of CPRNC binary passed in from CIME or specified in machine file. skip build
  find_program(CPRNC_BINARY cprnc ${CPRNC_DIR})
  IF ( CPRNC_BINARY )
     ADD_EXECUTABLE( cprnc IMPORTED)
     message("-- CPRNC_BINARY = ${CPRNC_BINARY}")
  ELSE()
     message(WARNING "cprnc not found in CPRNC_DIR")
     SET(CPRNC_DIR "")
  ENDIF ()
ENDIF ()

IF (NOT (CPRNC_DIR))
  # compile CPRNC from CIME source code. Requires CIME support for machine
  message("-- CPRNC_BINARY = will compile from source")
  message("-- If cmake aborts, set CPRNC_DIR to location of external cprnc executable")

  # cprnc's cmake wont search for netcdf, so we have to find it first:
  FIND_PACKAGE(NetCDF "4.0" COMPONENTS C Fortran)
  set(NETCDF_LIBRARIES ${NetCDF_Fortran_LIBRARIES})
  set(NETCDF_INCLUDE_DIR ${NetCDF_Fortran_INCLUDE_DIRS})

  # needed for CPRNC build system
  file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/utils/cime/tools/cprnc/Macros.cmake" "#dummy Macros file for non-CIME machines")

  SET (CPRNC_INSTALL_DIR ${CMAKE_BINARY_DIR}/utils/cime/tools/cprnc)
  SET (CPRNC_BINARY ${CMAKE_BINARY_DIR}/utils/cime/tools/cprnc/cprnc)
  ADD_SUBDIRECTORY(utils/cime/tools/cprnc)
ENDIF ()

IF (${HOMME_USE_KOKKOS})
  # Add unit tests for C++ code
  ADD_SUBDIRECTORY(test/unit_tests)
ENDIF ()

SET(EXEC_LIST ""  CACHE INTERNAL "List of configured executables")
IF(${BUILD_HOMME_SWEQX})
  ADD_SUBDIRECTORY(src/sweqx)
ENDIF()

IF(${BUILD_HOMME_PREQX})
  ADD_SUBDIRECTORY(src/preqx)
ENDIF()

IF(${BUILD_HOMME_THETA})
  ADD_SUBDIRECTORY(src/theta)
  ADD_SUBDIRECTORY(src/theta-l)
ENDIF()

IF(${BUILD_HOMME_PREQX_ACC})
  ADD_SUBDIRECTORY(src/preqx_acc)
ENDIF()

IF(${BUILD_HOMME_PREQX_KOKKOS})
  ADD_SUBDIRECTORY(src/preqx_kokkos)
ENDIF()

IF(${BUILD_HOMME_PESE})
  ADD_SUBDIRECTORY(src/pese)
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
