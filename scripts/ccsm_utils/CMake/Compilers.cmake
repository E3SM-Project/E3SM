# Flags for builds with different machines/compilers.
#
# This module is currently a catch-all for compiler-specific functionality
# needed by CESM. It defines OS and compiler CPP macros and CESM build
# types, as well as including the file containing CESM compiler flags, if
# necessary.
#
# There is also one function intended for CTest test writers, described
# below.
#
#==========================================================================
#
# define_Fortran_stop_failure
#
# Arguments:
#    test_name - Name of a CTest test.
#
# Ensures that if the named test uses "STOP 1" to signal failure, that this
# is detected by CTest. Currently this is only necessary for NAG, which
# prints the stop code rather than using it as an error code.
#
#==========================================================================

#==========================================================================
# Copyright (c) 2013-2014, University Corporation for Atmospheric Research
#
# This software is distributed under a two-clause BSD license, with no
# warranties, express or implied. See the accompanying LICENSE file for
# details.
#==========================================================================

#=================================================
# Add new build types.
#=================================================

# Add CESM build types.
set(CMAKE_Fortran_FLAGS_CESM "" CACHE STRING
  "Flags used by the Fortran compiler during CESM builds."
  FORCE)
set(CMAKE_C_FLAGS_CESM "" CACHE STRING
  "Flags used by the C compiler during CESM builds."
  FORCE)
mark_as_advanced(CMAKE_Fortran_FLAGS_CESM CMAKE_C_FLAGS_CESM)

set(CMAKE_Fortran_FLAGS_CESM_DEBUG "" CACHE STRING
  "Flags used by the Fortran compiler during CESM DEBUG builds."
  FORCE)
set(CMAKE_C_FLAGS_CESM_DEBUG "" CACHE STRING
  "Flags used by the C compiler during CESM DEBUG builds."
  FORCE)
mark_as_advanced(CMAKE_Fortran_FLAGS_CESM_DEBUG CMAKE_C_FLAGS_CESM_DEBUG)

set(all_build_types
  "None Debug Release RelWithDebInfo MinSizeRel CESM CESM_DEBUG")
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
  "Choose the type of build, options are: ${all_build_types}."
  FORCE)

#=================================================
# Define OS and compiler macros.
#=================================================

# Define OS.
string(TOUPPER ${CMAKE_SYSTEM_NAME} os)
add_definitions(-D${os})

# Define CESM-compatible compiler names.
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL NAG)
  set(compiler_name nag)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  set(compiler_name gnu)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL XL)
  set(compiler_name ibm)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  set(compiler_name intel)
endif()

# Define CPP macro for the compiler.
string(TOUPPER -DCPR${compiler_name} compiler_cppdef)
add_definitions(${compiler_cppdef})

#=================================================
# Utility functions.
#=================================================

# Add flags to space-separated list rather than normal CMake list.
function(add_flags list)
  string(REPLACE ";" " " flags "${ARGN}")
  set(${list} "${${list}} ${flags}" PARENT_SCOPE)
endfunction()

# Add configuration-specific preprocessor definitions.
function(add_config_definitions configuration)
  get_directory_property(cppdefs COMPILE_DEFINITIONS_${configuration})
  foreach(flag IN LISTS ARGN)
    string(REPLACE "-D" "" def "${flag}")
    list(APPEND cppdefs ${def})
  endforeach()
  set_directory_properties(PROPERTIES COMPILE_DEFINITIONS_${configuration}
    "${cppdefs}")
endfunction()

#=================================================
# Use CESM Macros file.
#=================================================

if("${CMAKE_BUILD_TYPE}" MATCHES CESM)
  include(${CMAKE_BINARY_DIR}/CESM_Macros.cmake)
endif()

#=================================================
# Build flags required to use pFUnit.
#=================================================

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  add_flags(CMAKE_Fortran_FLAGS -assume realloc_lhs)
endif()

#=================================================
# Add flags for debugging output.
#=================================================

# Define Fortran compiler flags.

# Add pretty output and extra warnings regardless of build type. However,
# don't set any options in the generic flags that would affect the
# generated binary, because we want to be able to get binaries that
# resemble what you get from the CESM flags.

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL NAG)
  add_flags(CMAKE_Fortran_FLAGS -strict95)
  if(USE_COLOR)
    add_flags(CMAKE_Fortran_FLAGS -colour)
  endif()

  # Add -kind=byte if it isn't anywhere else.
  if(NOT "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}" MATCHES -kind=byte)
    add_flags(CMAKE_Fortran_FLAGS -kind=byte)
  endif()

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  # Turn on warnings, but leave out uninitialized check as it was producing
  # a lot of false positives.
  add_flags(CMAKE_Fortran_FLAGS -Wall -Wextra -Wno-uninitialized)

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL XL)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
endif()

# Define C flags, analogous to the above Fortran block.
if(CMAKE_C_COMPILER_LOADED)
  if(${CMAKE_C_COMPILER_ID} STREQUAL GNU)
    add_flags(CMAKE_C_FLAGS -Wall -Wextra -pedantic)
  endif()
endif()

#=================================================
# Help CTest tests recognize when "stop X" is called with non-zero X.
#=================================================

# Detect "STOP" for CTest.
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL NAG)
  # NAG prints the stop code instead of yielding a non-zero return, so we
  # have to use a regex to catch that.
  function(define_Fortran_stop_failure test_name)
    set_tests_properties(${test_name} PROPERTIES
      FAIL_REGULAR_EXPRESSION "STOP: [1-9]")
  endfunction(define_Fortran_stop_failure)
else()
  # Usually, stop /= 0 is already detected with the return code.
  function(define_Fortran_stop_failure test_name)
  endfunction(define_Fortran_stop_failure)
endif()
