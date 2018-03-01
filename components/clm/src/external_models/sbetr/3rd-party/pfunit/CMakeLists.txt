# ------------------------------------------------------------------------ #
# DOCUMENTATION:
# ------------------------------------------------------------------------ #
#
# Command line options:
#     MPI=YES                                 ! defaults to NO
#     INSTALL_PATH=<path>                     ! defaults to ${PFUNIT}
#     MAX_RANK=<max array rank for generated code>
#                                             ! defaults to 5 or
#                                             !   ${PFUNIT_MAX_ARRAY_RANK}
# 
# Usage:
#   cmake -DMPI=YES -DINSTALL_PATH=<A path> <path to source>
#
# ------------------------------------------------------------------------ #
INCLUDE(CMakeForceCompiler)
cmake_minimum_required(VERSION 2.8)
project (pFUnit)
enable_language (Fortran)
 set (CMAKE_VERBOSE_MAKEFILE 1)

# MESSAGE(ENV{PFUNIT}:"$ENV{PFUNIT}") 
# set (ENV_PFUNIT "$ENV{PFUNIT}")
# MESSAGE(OK:"${ENV_PFUNIT}") 
# OPTION(INSTALL_PATH "Use PFUNIT for development" "$ENV{PFUNIT}")  
# OPTION(INSTALL_PATH "Use PFUNIT for development" "${ENV_PFUNIT}")  
OPTION(INSTALL_PATH "Use PFUNIT for development" OFF) 
IF(NOT INSTALL_PATH)
   IF(DEFINED ENV{PFUNIT})
      set (INSTALL_PATH "$ENV{PFUNIT}")
   ELSE()
      MESSAGE("** WARNING: Install path not found! Please set PFUNIT or use -DINSTALL_PATH=...")
   ENDIF()
ENDIF()
IF(NOT DEFINED ENV{PFUNIT})
   MESSAGE("++ Remember to set the environment variable PFUNIT to the install directory!")
   MESSAGE("++ Example: export PFUNIT=${INSTALL_PATH}")
ENDIF()
IF( NOT $ENV{PFUNIT} STREQUAL INSTALL_PATH )
   MESSAGE("** WARNING:  Environment variable PFUNIT not equal to INSTALL_PATH.")
   MESSAGE("**           Found: PFUNIT=$ENV{PFUNIT}")
ENDIF()
MESSAGE("-- Installing to INSTALL_PATH: ${INSTALL_PATH}")

OPTION(MAX_RANK "Maximum array rank for generated code." NO)
IF(NOT MAX_RANK)
   # Remove PFUNIT_MAX_RANK in version 4.
   IF(DEFINED ENV{PFUNIT_MAX_RANK})
      SET (MAX_RANK "$ENV{PFUNIT_MAX_RANK}")
   # Promote following to primary in version 4.
   ELSEIF(DEFINED ENV{PFUNIT_MAX_ARRAY_RANK})
      SET (MAX_RANK "$ENV{PFUNIT_MAX_ARRAY_RANK}")
   ELSE()
      SET (MAX_RANK 5)
   ENDIF()
ENDIF()


OPTION(MPI "Use MPI for parallel runs" NO)
OPTION(FIND_MPI "Attempt to automatically find MPI information" YES)
OPTION(MPI_USE_MPIEXEC "Use CMake's FindMPI to find mpiexec" NO)
OPTION(OPENMP "Use OPENMP for parallel runs" NO)

# The following sets robust to "ON".
OPTION(ROBUST "Build robust test runner" YES)

# If the user specifies this on the command line (or equivalently, in the
# CMake cache), we take them at their word and do not attempt any kind of
# autodetection.
set(OPENMP_FLAGS "OPENMP_FLAGS-NOTFOUND"
  CACHE STRING "OpenMP flags for the Fortran compiler.")


# OPTION doesn't seem to work as expected.
OPTION(INSTALL_PATH "Use PFUNIT for development" ENV{PFUNIT})
if ( NOT INSTALL_PATH )
   set (INSTALL_PATH "$ENV{PFUNIT}" )
endif()

if (MPI)
    add_definitions(-DUSE_MPI)
    if (FIND_MPI)
        find_package(MPI REQUIRED)
        include_directories(${MPI_INCLUDE_PATH})
        set (COMPILER_WRAPPER mpif90)
    endif()
    message( STATUS "MPI enabled")
endif()

# 2014-0807 MLR MAX_RANK needed for code generation, not here.
#if (MAX_RANK)
#   add_definitions(-DMAX_RANK=${MAX_RANK})
#endif()

if (ROBUST)
    add_definitions(-DBUILD_ROBUST)
endif()

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

    # Use ifort
    execute_process(COMMAND ifort --version 
                    OUTPUT_VARIABLE INTEL_COMPILER_VERSION)
    string (REGEX MATCH "([0-9]+)" 
            NUM_VERSION ${INTEL_COMPILER_VERSION})
    if("${NUM_VERSION}" STREQUAL "13")
        add_definitions(-DINTEL_13)
    endif()
    add_definitions(-DIntel)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g -check uninit -check bounds -traceback -assume realloc_lhs ${CMAKE_Fortran_FLAGS}")


elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)

    # Use Gfortran
    add_definitions(-DGNU)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fbounds-check")

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL NAG)

    # use nagfor			     
    add_definitions(-DNAG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -f2003 -w=uda -gline -fpp -mismatch_all -fpp -C=present")

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL XL)

    # Use XLF.
    add_definitions(-DIBM)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -WF,-qfpp -C")

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)

    # Use PGI.
    add_definitions(-DPGI)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g -traceback -Mallocatable=03 -Mbounds -Mchkfpstk -Mchkstk")

else()

    message( FATAL_ERROR "Unrecognized compiler. Please use ifort, gfortran, gfortran-mp-4.8, PGI, or nagfor" )

endif()

# Fortran OpenMP support is not yet integrated into the CMake distribution.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/include/cmake/Modules/")
if (OPENMP)
#    find_package(OpenMP_Fortran REQUIRED) 
#    if(OPENMP_FORTRAN_FOUND)
#      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#      set(CMAKE_Fortran_LINKER_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#      message( STATUS "OpenMP enabled")
#    endif()
   if (OPENMP_FLAGS)
     set(OpenMP_Fortran_FLAGS "${OPENMP_FLAGS}")
   else()
     find_package(OpenMP_Fortran REQUIRED)
   endif()
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
   set(CMAKE_Fortran_LINKER_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
   message( STATUS "OpenMP enabled")
endif()

set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
set(CMAKE_SKIP_RPATH ON)

add_subdirectory (source)
add_subdirectory (tests)
add_subdirectory (include)
add_subdirectory (bin)

# Override on command line with -DCMAKE_INSTALL_PREFIX=... or -DINSTALL_PATH=...
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${INSTALL_PATH}" CACHE PATH "default install path" FORCE )
endif()

# Packaging
set(CPACK_PACKAGE_CONTACT "pFUnit developers <pfunit-developers@list.sourceforge.net>")
include(CPack)
