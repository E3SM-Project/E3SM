## This is a CMake file, part of Unidata's netCDF Fortran package.
# Copyright 2012-2014, see the COPYRIGHT file for more information.

##################################
# Set Project Properties
##################################

#Minimum required CMake Version
cmake_minimum_required(VERSION 2.8.9)

#Project Name for netCDF-Fortran
PROJECT (NC4F Fortran C)
set(PACKAGE "netcdf-fortran" CACHE STRING "")
#Project Version
SET(NC4F_VERSION_MAJOR 4)
SET(NC4F_VERSION_MINOR 4)
SET(NC4F_VERSION_PATCH 3)
SET(NC4F_VERSION_NOTE "")
SET(NC4F_VERSION ${NC4F_VERSION_MAJOR}.${NC4F_VERSION_MINOR}.${NC4F_VERSION_PATCH}${NC4F_VERSION_NOTE})
SET(VERSION ${NC4F_VERSION})
SET(NC4F_LIB_VERSION 6.1.1)
SET(NC4F_SO_VERSION 6)
SET(PACKAGE_VERSION ${VERSION})

# Set known dependencies.
SET(EXTRA_DEPS "")
# Get system configuration,
# Use it to determine osname,
# os release, cpu. These will be used
# when committing to CDash.
find_program(UNAME NAMES uname)
IF(UNAME)
  macro(getuname name flag)
    exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
  endmacro(getuname)
  getuname(osname -s)
  getuname(osrel  -r)
  getuname(cpu    -m)
  set(TMP_BUILDNAME        "${osname}-${osrel}-${cpu}" CACHE STRING "Build name variable for CDash")
ENDIF()

###
# Allow for some customization of the buildname.
# This will make it easier to identify different builds,
# based on values passed from command line/shell scripts.
#
# For ctest scripts, we can use CTEST_BUILD_NAME.
###

SET(BUILDNAME_PREFIX "" CACHE STRING "")
SET(BUILDNAME_SUFFIX "" CACHE STRING "")

IF(BUILDNAME_PREFIX)
  SET(TMP_BUILDNAME "${BUILDNAME_PREFIX}-${TMP_BUILDNAME}")
ENDIF()

IF(BUILDNAME_SUFFIX)
  SET(TMP_BUILDNAME "${TMP_BUILDNAME}-${BUILDNAME_SUFFIX}")
ENDIF()

IF(NOT BUILDNAME)
  SET(BUILDNAME "${TMP_BUILDNAME}" CACHE STRING "Build name variable for CDash")
ENDIF()
###
# End BUILDNAME customization.
###


# Determine if DOXYGEN and DOT are available. These will be used
# when building the documentation.

FIND_PACKAGE(Doxygen)
FIND_PROGRAM(NC_DOT NAMES dot)


# For CMAKE_INSTALL_LIBDIR
INCLUDE(GNUInstallDirs)

IF(MSVC)
  SET(GLOBAL PROPERTY USE_FOLDERS ON)
ENDIF()

# auto-configure style checks, other CMake modules.
INCLUDE (${CMAKE_ROOT}/Modules/CheckLibraryExists.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckCXXSourceCompiles.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/TestBigEndian.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckSymbolExists.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/GetPrerequisites.cmake)
FIND_PACKAGE(PkgConfig QUIET)

# Enable 'dist and distcheck'.
# File adapted from http://ensc.de/cmake/FindMakeDist.cmake
# TODO: this can't currently be found, need
# FIND_PACKAGE(MakeDist)
# End 'enable dist and distcheck'

# Set the build type.
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING "Choose the type of build, options are: None, Debug, Release."
    FORCE)
ENDIF()

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

################################
# End Project Properties
################################


################################
# Set CTest Properties
################################

ENABLE_TESTING()
INCLUDE(CTest)

# Copy the CTest customization file into binary directory, as required.
# TODO: need to provide suitable CTestCustom.cmake


# Set Memory test program for non-MSVC based builds.
# Assume valgrind for now.
IF(NOT MSVC)
  SET (CTEST_MEMORYCHECK_COMMAND valgrind CACHE STRING "")
ENDIF()

# Set variable to define the build type.
INCLUDE(GenerateExportHeader)

################################
# End CTest Properties
################################


################################
# Compiler Configuration
################################


# Default to shared libs on.
OPTION(BUILD_SHARED_LIBS "Configure netCDF as a shared library." ON)
SET (LIB_TYPE STATIC)
IF (BUILD_SHARED_LIBS)
  SET(LIB_TYPE SHARED)
  IF(CMAKE_COMPILER_IS_GNUCC OR APPLE)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  ENDIF()
ENDIF()

# Supress unused variable and parameter warnings, for the time being,
# on GCC.
#
# Also, set some other default compiler flags.
IF(CMAKE_COMPILER_IS_GNUCC OR APPLE)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -Wno-unused-variable -Wno-unused-parameter")
ENDIF()

IF(NOT ENABLE_COVERAGE_TESTS)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
ENDIF()


ADD_DEFINITIONS()

# Supress CRT Warnings.
# Only necessary for Windows
IF(MSVC)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF()


################################
# End Compiler Configuration
################################


################################
# Define Utility Macros
################################

# Macro to append files to the EXTRA_DIST files.
SET(EXTRA_DIST "")
MACRO(ADD_EXTRA_DIST files)
  FOREACH(F ${files})
    SET(EXTRA_DIST ${EXTRA_DIST} ${CMAKE_CURRENT_SOURCE_DIR}/${F})
    SET(EXTRA_DIST ${EXTRA_DIST} PARENT_SCOPE)
  ENDFOREACH()
ENDMACRO()

# A basic script used to convert m4 files
FIND_PROGRAM(NC4F_M4 NAMES m4)
MACRO(GEN_m4 filename)
  IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c
      COMMAND ${NC4F_M4}
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c
      VERBATIM
      )
  ENDIF()
ENDMACRO(GEN_m4)

# [LEFT OFF HERE] replacing macros with ones from netcdf-c macros using NC4F_prefix

# Binary tests, but ones which depend on value of 'TEMP_LARGE' being defined.
MACRO(add_bin_env_temp_large_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()


  ADD_TEST(${prefix}_${F} bash "-c" "TEMP_LARGE=${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()


# Tests which are binary, but depend on a particular environmental variable.
MACRO(add_bin_env_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()


  ADD_TEST(${prefix}_${F} bash "-c" "TOPSRCDIR=${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()

# Build a binary used by a script, but don't make a test out of it.
MACRO(build_bin_test F)
  ADD_EXECUTABLE(${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()
ENDMACRO()

# Binary tests which are used by a script looking for a specific name.
MACRO(add_bin_test_no_prefix F)
  build_bin_test(${F})
  ADD_TEST(${F} ${EXECUTABLE_OUTPUT_PATH}/${F})
  IF(MSVC)
    SET_PROPERTY(TEST ${F} PROPERTY FOLDER "tests/")
  ENDIF()
ENDMACRO()

MACRO(add_bin_test prefix F)

  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${F}.F)
    set(ext ".F")
  ELSE()
    set(ext ".f90")
  ENDIF()

  ADD_EXECUTABLE(${prefix}_${F} ${F}${ext})
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()
  ADD_TEST(${prefix}_${F} ${EXECUTABLE_OUTPUT_PATH}/${prefix}_${F})
  IF(MSVC)
    SET_PROPERTY(TEST ${prefix}_${F} PROPERTY FOLDER "tests/")
  ENDIF()
ENDMACRO()

# A script to print out information at the end of the configuration step.
MACRO(print_conf_summary)
  MESSAGE("")
  MESSAGE("")
  MESSAGE(STATUS "CMake Summary:")
  MESSAGE("")
  MESSAGE(STATUS "\tBuild Type:\t\t\t" ${CMAKE_BUILD_TYPE})
  MESSAGE(STATUS "\tBuilding Shared Libraries:\t" ${BUILD_SHARED_LIBS})
  MESSAGE(STATUS "\tBuilding netCDF-4:\t\t" ${ENABLE_NETCDF4})
  MESSAGE(STATUS "\tBuilding DAP Support:\t\t" ${ENABLE_DAP})

  IF(STATUS_PNETCDF)
    MESSAGE(STATUS "\tUsing pnetcdf:\t\t\t" ${STATUS_PNETCDF})
  ENDIF(STATUS_PNETCDF)

  MESSAGE(STATUS "\tUsing Parallel IO:\t\t" ${STATUS_PARALLEL})

  IF(CMAKE_PREFIX_PATH)
    MESSAGE(STATUS "\tCMake Prefix Path:\t\t" "${CMAKE_PREFIX_PATH}")
  ENDIF()

  SET(ALL_DEP_LIBS "${NETCDF_C_LIBRARY};${EXTRA_DEPS}")
  MESSAGE(STATUS "\tLinking against:\t\t" "${ALL_DEP_LIBS}")
  FOREACH(_LIB ${EXTRA_DEPS})
    GET_FILENAME_COMPONENT(_LIB_NAME ${_LIB} NAME_WE)
    STRING(REGEX REPLACE "^lib" "" _NAME ${_LIB_NAME})
    LIST(APPEND NCF_LIBS "-l${_NAME}")
  ENDFOREACH()
  MESSAGE("")

ENDMACRO()
# Shell script Macro
MACRO(add_sh_test prefix F)
  IF(NOT MSVC)
    ADD_TEST(${prefix}_${F} bash "-c" "export srcdir=${CMAKE_CURRENT_SOURCE_DIR};export TOPSRCDIR=${CMAKE_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}/${F}.sh")
  ENDIF()

ENDMACRO()


# Macro for replacing '/MD' with '/MT'.
# Used only on Windows, /MD tells VS to use the shared
# CRT libs, MT tells VS to use the static CRT libs.
#
# Taken From:
# 	http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F
#
MACRO(specify_static_crt_flag)

  SET(vars
    CMAKE_C_FLAGS
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_MINSIZEREL
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELWITHDEBINFO)

  FOREACH(flag_var ${vars})
    IF(${flag_var} MATCHES "/MD")
      STRING(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    ENDIF()
  ENDFOREACH()

  FOREACH(flag_var ${vars})
    MESSAGE(STATUS " '${flag_var}': ${${flag_var}}")
  ENDFOREACH()
  MESSAGE(STATUS "")

ENDMACRO()

################################
# End Macro Definitions
################################


#Add custom CMake Module
SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/"
  CACHE INTERNAL "Location of our custom CMake modules.")

SET(USE_LOGGING "ON")

#####
# Some utility macros/scripts.

# Tests which are binary, but depend on a particular environmental variable.
MACRO(add_bin_env_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdff)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()


  ADD_TEST(${prefix}_${F} bash "-c" "TOPSRCDIR=${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()

# Build a binary used by a script, but don't make a test out of it.
MACRO(build_bin_test F ext)
  ADD_EXECUTABLE(${F} ${F}${ext})
  TARGET_LINK_LIBRARIES(${F} ncfortran netcdff ${NETCDF_C_LIBRARY})
  IF(MSVC)
    SET_TARGET_PROPERTIES(${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()
ENDMACRO()

OPTION(LARGE_FILE_TESTS "Run large file tests, which are slow and take lots of disk." OFF)
OPTION(BUILD_BENCHMARKS "Run F90 I/O Benchmarks" OFF)
OPTION(TEST_WITH_VALGRIND "Run extra tests with valgrind" OFF)
OPTION(TEST_PARALLEL "Run parallel I/O tests for F90 and F77" OFF)
IF(TEST_PARALLEL)
  SET(STATUS_PARALLEL ON)
ENDIF()

###
# Seek out dependent libraries.
###

IF(NOT netCDF_LIBRARIES AND NOT netCDF_INCLUDE_DIR)
  FIND_PACKAGE(netCDF QUIET)
ELSE()
  SET(netCDF_FOUND TRUE)
ENDIF()

IF (netCDF_FOUND)
  INCLUDE_DIRECTORIES(SYSTEM ${netCDF_INCLUDE_DIR})
  SET(NETCDF_C_LIBRARY ${netCDF_LIBRARIES})
  SET(NETCDF_C_INCLUDE_DIR ${netCDF_INCLUDE_DIR})
ELSE()
  # netCDF not properly packaged. Try to find it manually.
  FIND_LIBRARY(NETCDF_C_LIBRARY NAMES netcdf libnetcdf)
  IF(NOT NETCDF_C_LIBRARY)
    MESSAGE(FATAL_ERROR "libnetcdf not found. Please reinstall and try again.")
  ELSE()
    MESSAGE(STATUS "Found netcdf: ${NETCDF_C_LIBRARY}")
    FIND_PATH(NC_H_INCLUDE_DIR "netcdf.h")
    IF(NOT NC_H_INCLUDE_DIR)
      MESSAGE(FATAL_ERROR "Directory containing netcdf.h cannot be found. Please reinstall and try again.")
    ELSE()
      INCLUDE_DIRECTORIES(SYSTEM ${NC_H_INCLUDE_DIR})
    ENDIF()
  ENDIF()
ENDIF()

###
# End 'seek out dependent libraries'
###



###
# TODO: Replace these with automatically-determined
# checks, like in configure.ac.
###

CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nc_def_opaque "" USE_NETCDF4)
CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nccreate "" USE_NETCDF_V2)
CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nc_set_log_level "" USE_LOGGING)
CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} oc_open "" BUILD_DAP)
CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nc_use_parallel_enabled "" BUILD_PARALLEL)

##
# Set some status messages for the summary.
##

SET(STATUS_PARALLEL OFF)
IF(BUILD_PARALLEL)
  SET(STATUS_PARALLEL ON)
ENDIF()

IF(BUILD_DAP)
  SET(ENABLE_DAP ON)
ENDIF()

##
# End status messages/variables.
##

IF(USE_NETCDF4)
  SET(ENABLE_NETCDF4 ON)
  ADD_DEFINITIONS(-DUSE_NETCDF4)
  SET(ENABLE_NETCDF_4 ON CACHE BOOL "")
  SET(ENABLE_NETCDF4 ON CACHE BOOL "")
ENDIF()

IF (UNIX AND ${CMAKE_SIZEOF_VOID_P} MATCHES "8")
  ADD_DEFINITIONS(-DLONGLONG_IS_LONG)
  # for easier debugging of cfortran.h
  IF (APPLE)
  #  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmacro-backtrace-limit=0")
     SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftemplate-backtrace-limit=0")
  endif ()
ENDIF()

# Determine C/Fortran pointer compatibility.
try_compile(COMPILE_SUCCESS ${CMAKE_CURRENT_BINARY_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/CMakeExtras/test_c_ptrdiff_t.f90"
  )
IF (COMPILE_SUCCESS)
  ADD_DEFINITIONS(-DHAVE_TS29113_SUPPORT)
ELSE ()
  # If the compile fails, PTRDIFF_T is not defined.
  # Get sizeof of ptrdiff_t.
  SET(EXEC_NAME "${CMAKE_CURRENT_BINARY_DIR}/sizeof_ptrdiff_t")
  try_compile(COMPILE_SUCCESS2 ${CMAKE_CURRENT_BINARY_DIR}
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeExtras/sizeof_ptrdiff_t.c"
    COPY_FILE ${EXEC_NAME}
    )
  IF (COMPILE_SUCCESS2)
    EXECUTE_PROCESS(
      COMMAND ${EXEC_NAME}
      OUTPUT_VARIABLE SIZEOF_PTRDIFF_T
      RESULT_VARIABLE SIZEOF_RESULT
      )
    IF (SIZEOF_RESULT)
      MESSAGE (FATAL_ERROR "UNABLE TO DETERMINE SIZEOF PTRDIFF_T")
    ELSE (SIZEOF_RESULT)
      ADD_DEFINITIONS(-DSIZEOF_PTRDIFF_T=${SIZEOF_PTRDIFF_T})
    ENDIF (SIZEOF_RESULT)
  ELSE()
    MESSAGE(FATAL_ERROR "Unable to compile ptrdiff")
  ENDIF()
endif()

IF(USE_NETCDF_V2)
  SET(BUILD_V2 ON CACHE BOOL "")
ENDIF()
# Turn this on by default when it's working
OPTION(ENABLE_FORTRAN_TYPE_CHECKS
  "Determine Fortran types corresponding to netCDF types" OFF)

# Toggle whether or not to run tests.
OPTION(ENABLE_TESTS "Enable netcdf-fortran tests." ON)

# Support old version-2 API for backward compatibility
SET (ENABLE_LOGGING "ON")
IF(ENABLE_LOGGING)
  SET(USE_LOGGING ON CACHE BOOL "")
ENDIF()

# Set the default fortran builds; default is to build f03
SET(BUILD_F90 "ON")
SET(BUILD_V2 "ON")
SET(BUILD_F03 "OFF")

IF(ENABLE_FORTRAN_TYPE_CHECKS)
  # Do tests to determine which Fortran types correspond to NCBYTE, NCSHORT, ...
  # The configure file got these by testing an F77 program, invoking
  # UD_FORTRAN_TYPES defined in acinclude.m4.  TODO: check the types
  # from the test program instead of assuming these defaults.
ELSE()
  # Set Fortran types to default.
  SET(NCBYTE_T "byte")
  SET(NCBYTE_T "byte")
  SET(NCSHORT_T "integer*2")
  SET(NF_INT1_T "byte")
  SET(NF_INT2_T "integer*2")
  # Default is for following to be true.  TODO: test instead
  SET(NF_INT1_IS_C_SIGNED_CHAR "1")
  SET(NF_INT2_IS_C_SHORT "1")
  SET(NF_INT_IS_C_INT "1")
  SET(NF_REAL_IS_C_FLOAT "1")
  SET(NF_DOUBLEPRECISION_IS_C_DOUBLE "1")
ENDIF(ENABLE_FORTRAN_TYPE_CHECKS)

# Create libsrc/nfconfig.inc (don't need nfconfig1.inc temporary file used with autoconf)
configure_file("${NC4F_SOURCE_DIR}/libsrc/nfconfig.in.cmake"
  "${NC4F_BINARY_DIR}/libsrc/nfconfig.inc")
# Create config.h file
configure_file("${NC4F_SOURCE_DIR}/config.h.cmake.in"
  "${NC4F_BINARY_DIR}/config.h")

# For now, just copy a stub file
FILE(COPY "${NC4F_SOURCE_DIR}/nf-config.cmake.in"
  DESTINATION "${NC4F_BINARY_DIR}"
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
FILE(RENAME "${NC4F_BINARY_DIR}/nf-config.cmake.in" "${NC4F_BINARY_DIR}/nf-config")

# make sure previous two files get cleaned up...
SET_DIRECTORY_PROPERTIES (DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${netcdf-fortran_BINARY_DIR}/config.h ${netcdf-fortran_BINARY_DIR}/nf-config)

##
# Configuration for post-install RPath
# Adapted from http://www.cmake.org/Wiki/CMake_RPATH_handling
##
IF(NOT MSVC)
  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  IF(APPLE)
    SET(CMAKE_MACOSX_RPATH ON)
  ENDIF()

  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  ENDIF("${isSystemDir}" STREQUAL "-1")
ENDIF()

##
# End configuration for post-install RPath
##

###
# Test-related options
###

# Set CDash-related options.
SET(NC_CTEST_PROJECT_NAME "netcdf-fortran" CACHE STRING "Project Name for CTest-based testing purposes.")
SET(NC_CTEST_DROP_SITE "my.cdash.org" CACHE STRING "Dashboard location for CTest-based testing purposes.")
SET(NC_CTEST_DROP_LOC_PREFIX "" CACHE STRING "Prefix for Dashboard location on remote server when using CTest-based testing.")

# Generate CTestConfig file from template.
CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/CTestConfig.cmake.in"
  "${CMAKE_SOURCE_DIR}/CTestConfig.cmake"
  @ONLY)


###
# Allow the user to specify libraries
# to link against, similar to automakes 'LIBS' variable.
###
SET(NC_EXTRA_DEPS "" CACHE STRING "Additional libraries to link against.")
IF(NC_EXTRA_DEPS)
  STRING(REPLACE " " ";" DEPS_LIST ${NC_EXTRA_DEPS})
  FOREACH(_DEP ${DEPS_LIST})
    STRING(REGEX REPLACE "^-l" "" _LIB ${_DEP})
    FIND_LIBRARY("${_LIB}_DEP" NAMES "${_LIB}" "lib${_LIB}")
    MESSAGE(${${_LIB}_DEP})
    IF(NOT "${_LIB}_DEP")
      MESSAGE(FATAL_ERROR "Error finding ${_LIB}.")
    ELSE()
      MESSAGE(STATUS "Found ${_LIB}: ${${_LIB}_DEP}")
    ENDIF()
    SET(EXTRA_DEPS ${EXTRA_DEPS} "${${_LIB}_DEP}")
  ENDFOREACH()
  MESSAGE("Extra deps: ${EXTRA_DEPS}")
  LIST(REMOVE_DUPLICATES EXTRA_DEPS)
ENDIF()
###
# End user-specified dependent libraries.
###

# Determine whether or not to generate documentation.
OPTION(ENABLE_DOXYGEN "Enable generation of doxygen-based documentation." OFF)
IF(ENABLE_DOXYGEN)
  IF(DOXYGEN_FOUND)
    # Offer the option to build internal documentation.
    OPTION(ENABLE_INTERNAL_DOCS "Build internal documentation. This is of interest to developers only." OFF)
    IF(ENABLE_INTERNAL_DOCS)
      SET(BUILD_INTERNAL_DOCS YES CACHE STRING "")
    ELSE()
      SET(BUILD_INTERNAL_DOCS NO CACHE STRING "")
    ENDIF()

    # Option to turn on the TODO list in the doxygen-generated documentation.
    OPTION(ENABLE_DOXYGEN_TASKS "Turn on test, todo, bug lists in documentation. This is of interest to developers only." OFF)
    IF(ENABLE_DOXYGEN_TASKS)
      SET(SHOW_DOXYGEN_TAG_LIST YES CACHE STRING "")
    ELSE(ENABLE_DOXYGEN_TASKS)
      SET(SHOW_DOXYGEN_TODO_LIST NO CACHE STRING "")
    ENDIF(ENABLE_DOXYGEN_TASKS)

    # Specify whether or not 'dot' was found on the system path.
    IF(NC_DOT)
      SET(HAVE_DOT YES CACHE STRING "")
    ELSE(NC_DOT)
      SET(HAVE_DOT NO CACHE_STRING "")
    ENDIF(NC_DOT)

  ELSE(DOXYGEN_FOUND)
    MESSAGE(STATUS "Unable to build internal documentation.  Doxygen does not appear to be on your executable path. Install doxygen and configure the project again.")
    SET(ENABLE_DOXYGEN OFF)
  ENDIF(DOXYGEN_FOUND)
ENDIF()

###
# Set global linker flags
###
IF(CMAKE_LINK_FLAGS)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_LINK_FLAGS}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_LINK_FLAGS}")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${CMAKE_LINK_FLAGS}")
ENDIF()


INSTALL(PROGRAMS ${NC4F_BINARY_DIR}/nf-config
        DESTINATION "${CMAKE_INSTALL_BINDIR}"
        COMPONENT utilities)

INCLUDE_DIRECTORIES(${NC4F_BINARY_DIR})
# End autotools-style checks for config.h

# Recurse into other subdirectories.
ADD_SUBDIRECTORY(fortran)

ADD_SUBDIRECTORY(libsrc)

IF(ENABLE_TESTS)
  ADD_SUBDIRECTORY(nf_test)
  ADD_SUBDIRECTORY(nf03_test)
ENDIF()

# Build the Fortran examples
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(docs)
print_conf_summary()
