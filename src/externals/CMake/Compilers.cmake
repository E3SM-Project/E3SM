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
