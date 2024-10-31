IF (NOT DEFINED PROJECT_NAME)
  SET (HOMME_STANDALONE TRUE CACHE BOOL "Whether Homme is built in standalone mode" FORCE)

  # Best practice to put cmake required version at the top
  CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

  # Print the sha of the last commit,
  # useful to double check which version was build in the build folder.
  EXECUTE_PROCESS (COMMAND git rev-parse HEAD
                     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                     OUTPUT_VARIABLE LAST_GIT_COMMIT_SHA
                     OUTPUT_STRIP_TRAILING_WHITESPACE)
  SET (LAST_GIT_COMMIT_SHA ${LAST_GIT_COMMIT_SHA})
  MESSAGE (STATUS "The sha of the last commit is ${LAST_GIT_COMMIT_SHA}")

ELSE ()
  SET (HOMME_STANDALONE FALSE CACHE BOOL "Whether Homme is built in standalone mode" FORCE)

  # Manually enable languages, in case they were not enabled by host project
  ENABLE_LANGUAGE(C)
  ENABLE_LANGUAGE(Fortran)
ENDIF ()

# By default, build executables only in standalone mode
OPTION (HOMME_BUILD_EXECS "Whether Homme should build executables and tests." ${HOMME_STANDALONE})

PROJECT(HOMME C Fortran)

INCLUDE(FortranCInterface)

IF (${HOMME_BUILD_EXECS})
  OPTION(HOMME_ENABLE_TESTING "Whether tests should be enabled" TRUE)
  IF (${HOMME_ENABLE_TESTING})
    # CTest For Nightly/Experimental Builds
    OPTION(HOMME_USE_CTEST "Use CTest Dashboard" FALSE)
    IF (${HOMME_USE_CTEST})
      INCLUDE(CTest)
    ENDIF ()

    # For testing
    ENABLE_TESTING()
    SET (HOMME_TESTING_PROFILE "nightly" CACHE STRING "Determine how long/pervasive the testing is. Currently available options: 'dev', 'short', 'nightly'.")
  ENDIF()
ENDIF()

# Option to build without pio and any netcdf in homme (no mesh runs either)
# The default is not on (meaning HAVE pio support) if we build any exec.
IF (HOMME_BUILD_EXECS OR BUILD_HOMME_TOOL)
  SET(BUILD_HOMME_WITHOUT_PIOLIBRARY OFF CACHE BOOL "")
ELSE ()
  SET(BUILD_HOMME_WITHOUT_PIOLIBRARY ON CACHE BOOL "")
ENDIF ()

SET(HOMME_BUILD_SCORPIO ON)
IF(BUILD_HOMME_WITHOUT_PIOLIBRARY)
  MESSAGE(STATUS "This configuration builds without PIO library and NetCDF calls")
  SET(HOMME_BUILD_SCORPIO OFF)
ELSE()
  IF (TARGET pioc OR TARGET piof)
    message(STATUS "pioc and/or piof are already defined; skipping SCORPIO build in HOMME")
    SET(HOMME_BUILD_SCORPIO OFF)
  ENDIF()
ENDIF()

#option to discard forcings call and push to/from F
#set it on in cmake file
SET(BUILD_HOMMEXX_BENCHMARK_NOFORCING OFF CACHE BOOL "")
if(BUILD_HOMMEXX_BENCHMARK_NOFORCING)
  message(STATUS "Build without test forcing and pushes to/from F, PERFORMANCE ONLY")
else()
  message(STATUS "Option BUILD_HOMMEXX_BENCHMARK_NOFORCING is OFF")
endif()

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
#  ${HOMME_BINARY_DIR}/tests/compilerTests/
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
OPTION(BUILD_HOMME_THETA_KOKKOS "SB81 NH model with theta and Exner pressure, Kokkos target" OFF)
OPTION(BUILD_HOMME_PREQX_ACC "Primitive equations FEM with OpenACC" ON)
OPTION(BUILD_HOMME_PREQX_KOKKOS "Primitive equations FEM with Kokkos" OFF)
OPTION(BUILD_HOMME_SWIM  "Shallow water equations implicit" OFF)
OPTION(BUILD_HOMME_PRIM  "Primitive equations implicit" OFF)
OPTION(BUILD_HOMME_TOOL  "Offline tool" ON)
OPTION(HOMME_ENABLE_COMPOSE "Build COMPOSE semi-Lagrangian tracer transport code" ON)
OPTION(HOMME_USE_SCORPIO  "Use Scorpio as the I/O library (Disable to use Scorpio classic)" ON)

#by default we don't need cxx
SET(HOMME_USE_CXX FALSE)
SET(HOMME_USE_KOKKOS FALSE)

# swim and prim executables require C++ compiler
IF(${BUILD_HOMME_SWIM} OR ${BUILD_HOMME_PRIM})
  SET (HOMME_USE_CXX TRUE)
  SET (HOMME_USE_TRILINOS TRUE)
ENDIF ()
#if trilinos is needed, we need cxx
IF (HOMME_USE_TRILINOS)
  SET (HOMME_USE_CXX TRUE)
ENDIF()
#if kokkos is needed, then so is cxx
SET (HOMME_USE_KOKKOS FALSE)
SET (HOMME_BUILD_KOKKOS_DYCORE FALSE)

IF (BUILD_HOMME_PREQX_KOKKOS OR BUILD_HOMME_THETA_KOKKOS)
  SET (HOMME_BUILD_KOKKOS_DYCORE TRUE)
ENDIF ()
IF (HOMME_BUILD_KOKKOS_DYCORE)
  SET (HOMME_USE_KOKKOS TRUE)

  # This option will only help propagate tests, it won't make bfb builds.
  # Compiler flags should be set separately in *.cmake file.
  IF (HOMME_ENABLE_TESTING)
    OPTION (HOMMEXX_BFB_TESTING "Whether we want bfb comparison cpu-vs-gpu and f90-vs-cxx" OFF)
  ENDIF ()

  SET (HOMME_USE_CXX TRUE)
  IF (NOT BUILD_HOMME_PREQX AND BUILD_HOMME_PREQX_KOKKOS AND HOMMEXX_BFB_TESTING)
    # If we build preqx kokkos, we also build preqx, so we can compare
    MESSAGE (STATUS "Setting manually disabled BUILD_HOMME_PREQX to ON, since BUILD_HOMME_PREQX_KOKKOS is ON, "
                    "and HOMMEXX_BFB_TESTING=ON.")

    # We need to compare against F90 implementation, so turn this on
    SET (BUILD_HOMME_PREQX TRUE)
  ENDIF ()

  IF (NOT BUILD_HOMME_THETA AND BUILD_HOMME_THETA_KOKKOS AND HOMMEXX_BFB_TESTING)
    # If we build theta kokkos, we also build theta, so we can compare
    MESSAGE (STATUS "Setting manually disabled BUILD_HOMME_THETA to ON, since BUILD_HOMME_THETA_KOKKOS is ON, "
                    "and HOMMEXX_BFB_TESTING=ON.")

    # We need to compare against F90 implementation, so turn this on
    SET (BUILD_HOMME_THETA TRUE)
  ENDIF ()
ENDIF ()

IF (BUILD_HOMME_PREQX OR BUILD_HOMME_THETA OR BUILD_HOMME_PREQX_ACC)
  SET (HOMME_BUILD_FORTRAN_DYCORE TRUE)
ENDIF()

IF (HOMME_ENABLE_COMPOSE)
  SET (HOMME_USE_CXX TRUE)
  SET (HOMME_USE_KOKKOS TRUE)
ENDIF ()

IF (HOMME_USE_CXX)
  MESSAGE(STATUS "This configuration of HOMME requires a C++ compiler")
  ENABLE_LANGUAGE(CXX)
ENDIF ()

# Option to use ARKode package from SUNDIALS
OPTION(HOMME_USE_ARKODE "Use ARKode package from SUNDIALS" FALSE)
IF (${HOMME_USE_ARKODE})
  MESSAGE(STATUS "Using ARKode, adding -DARKODE")
  ADD_DEFINITIONS(-DARKODE)
ENDIF ()

IF (HOMME_USE_KOKKOS)
  ############################################
  # Selection of Kokkos execution space
  ############################################

  SET (HOMMEXX_EXEC_SPACE "Default" CACHE STRING "Select the kokkos exec space")

  STRING (TOUPPER ${HOMMEXX_EXEC_SPACE} HOMMEXX_EXEC_SPACE_UPPER)

  IF (${HOMMEXX_EXEC_SPACE_UPPER} STREQUAL "SYCL")
    SET (HOMMEXX_SYCL_SPACE ON)
  ELSEIF (${HOMMEXX_EXEC_SPACE_UPPER} STREQUAL "HIP")
    SET (HOMMEXX_HIP_SPACE ON)
  ELSEIF (HOMMEXX_EXEC_SPACE_UPPER STREQUAL "CUDA")
    SET (HOMMEXX_CUDA_SPACE ON)
  ELSEIF (HOMMEXX_EXEC_SPACE_UPPER STREQUAL "OPENMP")
    SET (HOMMEXX_OPENMP_SPACE ON)
  ELSEIF (HOMMEXX_EXEC_SPACE_UPPER STREQUAL "THREADS")
    SET (HOMMEXX_THREADS_SPACE ON)
  ELSEIF (HOMMEXX_EXEC_SPACE_UPPER STREQUAL "SERIAL")
    SET (HOMMEXX_SERIAL_SPACE ON)
  ELSEIF (HOMMEXX_EXEC_SPACE_UPPER STREQUAL "DEFAULT")
    SET (HOMMEXX_DEFAULT_SPACE ON)
  ELSE()
    MESSAGE (ABORT "Invalid choice for 'HOMMEXX_EXEC_SPACE'. Valid options (case insensitive) are 'Cuda', 'OpenMP', 'Threads', 'Serial', 'Default'")
  ENDIF()

  # Execution space parameters. 8 is a nice size for V100.
  SET (HOMMEXX_CUDA_MIN_WARP_PER_TEAM 8 CACHE STRING "Minimum number of warps to get 100% occoupancy on GPU")
  SET (HOMMEXX_CUDA_MAX_WARP_PER_TEAM 16 CACHE STRING "Maximum number of warps to get 100% occoupancy on GPU")

  # An option to allow to use GPU pointers for MPI calls. The value of this option is irrelevant for CPU/KNL builds.
  OPTION (HOMMEXX_MPI_ON_DEVICE "Whether we want to use device pointers for MPI calls (relevant only for GPU builds)" ON)

  # An option to allow workspace sharing on GPU
  OPTION (HOMMEXX_CUDA_SHARE_BUFFER "Whether we want to allow for buffer sharing on GPU. This feature incurs some computational overhead but can allow running of larger problems (relevant only for GPU builds)" OFF)
ENDIF()

##############################################################################
# Compiler options set in cmake/SetCompiler.cmake
##############################################################################
IF (HOMME_STANDALONE)
  INCLUDE(${HOMME_SOURCE_DIR}/cmake/SetCompilerFlags.cmake)
ENDIF()

# Create CsmShare target needed for test execs. If Homme is built
# inside another project, then the host project is in charge of
# ensuring that a csm_share target exists at this point
if (HOMME_STANDALONE)
  add_subdirectory(utils/csm_share)
elseif (NOT TARGET csm_share)
  string (CONCAT msg
    "Error! Homme needs to link against csm_share libraries.\n"
    "       Please ensure that a csm_share target is built by\n"
    "       the host project *before* entering the homme directory.\n")
  message ("${msg}")
  message (FATAL_ERROR "Aborting.")
endif()

OPTION(HOMME_USE_MKL "Whether to use Intel's MKL/oneMKL instead of blas/lapack" FALSE)
IF(HOMME_USE_MKL)
  IF(DEFINED ENV{MKLROOT})
    SET(MKL_ROOT_PATH "$ENV{MKLROOT}")
    IF(MKL_ROOT_PATH MATCHES "oneapi/mkl")
      SET(MKL_TYPE "oneMKL")
      MESSAGE(STATUS "Detected oneMKL based on MKLROOT: ${MKL_ROOT_PATH}")
      FIND_PACKAGE(MKL REQUIRED $ENV{MKLROOT}/lib/cmake/mkl)
    ELSE()
      SET(MKL_TYPE "Intel MKL")
      MESSAGE(STATUS "Detected standalone Intel MKL based on MKLROOT: ${MKL_ROOT_PATH}")
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR "MKLROOT environment variable is not set. Please set it to your MKL installation path.")
  ENDIF()
  MESSAGE(STATUS "HOMME_USE_MKL is ON. The flag -mkl/-qmkl will be added to each executable/library.")
ELSE()
  OPTION(HOMME_FIND_BLASLAPACK "Whether to use system blas/lapack" FALSE)
  MESSAGE(STATUS "HOMME_FIND_BLASLAPACK=${HOMME_FIND_BLASLAPACK}")

  IF(HOMME_FIND_BLASLAPACK)
    FIND_PACKAGE(BLAS REQUIRED)
    FIND_PACKAGE(LAPACK REQUIRED)

    # Set libraries in cache, so they are accessible also from outside the Homme project.
    # That's the case if the host project adds the homme subdir from folder foo/bar1, but
    # then calls createExecLib from foo/bar2 (non cache variables would not survive outside bar1)
    # Note: an alternative would be to wrap them in a cmake target, but that's more complicated.
    SET (BLAS_LIBRARIES "${BLAS_LIBRARIES}" CACHE INTERNAL "")
    SET (LAPACK_LIBRARIES "${LAPACK_LIBRARIES}" CACHE INTERNAL "")
  ELSE()
    SET(HOMME_FIND_BLASLAPACK FALSE)
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

SET (HOMMEXX_ENABLE_GPU FALSE)
SET (HOMMEXX_ENABLE_GPU_F90 FALSE) 

IF (HOMME_USE_KOKKOS)

  IF (CUDA_BUILD OR HIP_BUILD OR SYCL_BUILD)
    SET (DEFAULT_VECTOR_SIZE 1)
    SET (HOMMEXX_ENABLE_GPU TRUE)
    SET (HOMMEXX_ENABLE_GPU_F90 TRUE) 
  ELSE ()
    SET (DEFAULT_VECTOR_SIZE 8)
  ENDIF()

  SET (HOMMEXX_VECTOR_SIZE ${DEFAULT_VECTOR_SIZE} CACHE STRING
	  "If AVX or Cuda or HIP or SYCL don't take priority, use this software vector size.")

  IF (CMAKE_BUILD_TYPE_UPPER MATCHES "DEBUG" OR CMAKE_BUILD_TYPE_UPPER MATCHES "RELWITHDEBINFO")
    SET (HOMMEXX_DEBUG ON)
  ENDIF()

  CONFIGURE_FILE (${CMAKE_CURRENT_SOURCE_DIR}/src/share/cxx/Hommexx_config.h.in
                  ${HOMME_BINARY_DIR}/src/share/cxx/Hommexx_config.h)

ENDIF ()

CONFIGURE_FILE (${CMAKE_CURRENT_SOURCE_DIR}/cmake/homme_git_sha.h.in ${HOMME_BINARY_DIR}//homme_git_sha.h)


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

# this section is for linking to SUNDIALS/ARKODE as an external library
IF (DEFINED SUNDIALS_DIR)
  LINK_DIRECTORIES(${SUNDIALS_DIR}/lib64)
  INCLUDE_DIRECTORIES(${SUNDIALS_DIR}/include)
  INCLUDE_DIRECTORIES(${SUNDIALS_DIR}/fortran)
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
ADD_SUBDIRECTORY(utils/cime/CIME/non_py/src/timing)


# CMAKE_CURRENT_SOURCE_DIR == homme
SET (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
SET (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/utils/cime/src/CMake" ${CMAKE_MODULE_PATH})
if (HOMME_BUILD_SCORPIO)
  SET (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/utils/externals/scorpio/cmake" ${CMAKE_MODULE_PATH})
endif ()

#
# Scorpio
#
#set search paths for SCORPIO's findNetCDF
# HOMME machine files set NETCDF_DIR, PNETCDF_DIR which we copy
# to variables used by SCORPIO.  Newer machine files sould direclty set
# necessary SCORPIO variables:
#
ADD_DEFINITIONS(-D_NO_MPI_RSEND)

if (HOMME_BUILD_SCORPIO)
  SET(NetCDF_PATH ${NETCDF_DIR})
  SET(PnetCDF_PATH ${PNETCDF_DIR})
  message("-- NetCDF_PATH = ${NetCDF_PATH}")
  message("-- PnetCDF_PATH = ${PnetCDF_PATH}")
  message("-- NetCDF_Fortran_PATH = ${NetCDF_Fortran_PATH}")
  message("-- NetCDF_C_PATH = ${NetCDF_C_PATH}")
  message("-- PnetCDF_Fortran_PATH = ${PnetCDF_Fortran_PATH}")
  message("-- PnetCDF_C_PATH = ${PnetCDF_C_PATH}")
  message("** Configuring SCORPIO")
  # pio needs cime/externals/genf90/genf90.pl
  if (HOMME_USE_SCORPIO)
    # Need to use scorpio's older genf90 for now
    SET(GENF90_PATH ${CMAKE_CURRENT_SOURCE_DIR}/utils/utils/externals/scorpio/src/genf90)
    SET(PIO_ENABLE_TOOLS OFF CACHE BOOL "Disabling Scorpio tool build")
    ADD_SUBDIRECTORY(utils/externals/scorpio)
    SET(PIO_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/utils/externals/scorpio/src/flib")
  else ()
    # The default I/O library used in "Scorpio classic"
    ADD_SUBDIRECTORY(utils/externals/scorpio_classic)
    SET(GENF90_PATH ${CMAKE_CURRENT_SOURCE_DIR}/utils/cime/CIME/non_py/externals/genf90)
    SET(PIO_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/utils/externals/scorpio/pio")
  endif ()
  INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS})
endif ()


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

#compile cprnc only if pio is built
if(HOMME_BUILD_EXECS AND NOT BUILD_HOMME_WITHOUT_PIOLIBRARY)
  IF (NOT (CPRNC_DIR))
    # compile CPRNC from CIME source code. Requires CIME support for machine
    message("-- CPRNC_BINARY = will compile from source")
    message("-- If cmake aborts, set CPRNC_DIR to location of external cprnc executable")

    # cprnc's cmake wont search for netcdf, so we have to find it first:
    FIND_PACKAGE(NetCDF "4.0" COMPONENTS C Fortran)
    set(NETCDF_LIBRARIES ${NetCDF_Fortran_LIBRARIES})
    set(NETCDF_INCLUDE_DIR ${NetCDF_Fortran_INCLUDE_DIRS})

    # needed for CPRNC build system
    set (cprnc_dummy_file "${CMAKE_CURRENT_SOURCE_DIR}/utils/cime/CIME/non_py/cprnc/Macros.cmake")
    if (NOT EXISTS "${cprnc_dummy_file}")
      file(WRITE "${cprnc_dummy_file}" "#dummy Macros file for non-CIME machines")
    endif ()

    # cprnc's CMake system needs various nonstandard variables:
    set(FFLAGS ${CMAKE_Fortran_FLAGS})
    set(BLDROOT ${CMAKE_CURRENT_SOURCE_DIR}/utils/cime/CIME/non_py/cprnc)
    set(SFC ${CMAKE_Fortran_COMPILER})
    set(SCC ${CMAKE_C_COMPILER})
    
    SET (CPRNC_INSTALL_DIR ${HOMME_BINARY_DIR}/utils/cime/CIME/non_py/cprnc)
    SET (CPRNC_BINARY ${HOMME_BINARY_DIR}/utils/cime/CIME/non_py/cprnc/cprnc)
    ADD_SUBDIRECTORY(utils/cime/CIME/non_py/cprnc)
  ENDIF ()
ENDIF ()

IF (HOMME_USE_KOKKOS)
  # Add ekat's cmake scripts folders to cmake path
  set (EKAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../externals/ekat)
  set (EKAT_CMAKE_PATH ${EKAT_SOURCE_DIR}/cmake)
  list(APPEND CMAKE_MODULE_PATH
    ${EKAT_CMAKE_PATH}
    ${EKAT_CMAKE_PATH}/tpls
  )

  # We first try to use find_package. If that doesn't work, build from EKAT's submodule
  include (EkatFindKokkos)
  if (NOT Kokkos_FOUND)
    # The following script checks if Kokkos is already available as a target, and if so does nothing.
    # For instance, if HOMME is built inside EAMxx, Kokkos will already be available
    include (EkatBuildKokkos)
  endif()
ENDIF ()

# This folder contains the CMake macro used to build cxx unit tests
# Add unit tests for C++ code
ADD_SUBDIRECTORY(test/unit_tests)

# Add libcompose, consisting of just C++ files, so we can build it just once for
# all exes. It does not depend on compile-time sizes.
if (HOMME_ENABLE_COMPOSE)
  message (STATUS "COMPOSE semi-Lagrangian transport enabled")
  if (HOMME_BUILD_KOKKOS_DYCORE)
    set (HOMME_BUILD_COMPOSE_CPP TRUE)
  endif()
  if (HOMME_BUILD_FORTRAN_DYCORE OR BUILD_HOMME_TOOL)
    # Fortran dycore(s) and the homme_tool exec do need composef90
    set (HOMME_BUILD_COMPOSE_F90 TRUE)
  endif()
else ()
  message (STATUS "COMPOSE semi-Lagrangian transport was explicitly disabled")
endif()

if (HOMME_BUILD_COMPOSE_F90)
  set (COMPOSE_LIBRARY_F90 "composef90")
  set (COMPOSE_LIBRARY ${COMPOSE_LIBRARY_F90})
  set (COMPOSE_PORT FALSE)
  add_subdirectory(src/share/compose ${COMPOSE_LIBRARY_F90})
endif()

if (HOMME_BUILD_COMPOSE_CPP)
  set (COMPOSE_LIBRARY_CPP "composec++")
  set (COMPOSE_LIBRARY ${COMPOSE_LIBRARY_CPP})
  set (COMPOSE_PORT TRUE)
  add_subdirectory(src/share/compose ${COMPOSE_LIBRARY_CPP})
endif ()

SET(EXEC_LIST ""  CACHE INTERNAL "List of configured executables")
IF(${BUILD_HOMME_SWEQX})
  ADD_SUBDIRECTORY(src/sweqx)
ENDIF()

IF(${BUILD_HOMME_PREQX})
  ADD_SUBDIRECTORY(src/preqx)
ENDIF()

IF(${BUILD_HOMME_THETA})
  ADD_SUBDIRECTORY(src/theta-l)
ENDIF()

IF(${BUILD_HOMME_PREQX_ACC})
  ADD_SUBDIRECTORY(src/preqx_acc)
ENDIF()

IF(${BUILD_HOMME_PREQX_KOKKOS})
  ADD_SUBDIRECTORY(src/preqx_kokkos)
ENDIF()

IF(${BUILD_HOMME_THETA_KOKKOS})
  ADD_SUBDIRECTORY(src/theta-l_kokkos)
ENDIF()
IF(${BUILD_HOMME_SWIM})
  ADD_SUBDIRECTORY(src/swim)
ENDIF()

IF(${BUILD_HOMME_PRIM})
  ADD_SUBDIRECTORY(src/prim)
ENDIF()

IF(${BUILD_HOMME_TOOL})
  ADD_SUBDIRECTORY(src/tool)
ENDIF()

# Tests and execs only if this is a standalone build
IF (${HOMME_BUILD_EXECS})
  # Set a default value for USE_MPIEXEC.
  # NOTE: if USE_MPIEXEC is already set, this line has no effect.
  SET (USE_MPIEXEC "mpiexec" CACHE STRING "")

  ADD_SUBDIRECTORY(test_execs)
  ADD_SUBDIRECTORY(dcmip_tests)
ENDIF()

# Create a "make clean-all" target to delete CMake cache and the cprnc build
ADD_CUSTOM_TARGET(clean-all
   COMMAND ${CMAKE_BUILD_TOOL} clean
   COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/clean-cmake-files.cmake"
)
