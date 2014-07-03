# Module used for CESM testing.
#
# This module contains statements that would otherwise be boilerplate in
# most CESM tests. It enables CTest testing, handles the USE_COLOR and
# ENABLE_GENF90 arguments, and includes several other modules.

#==========================================================================
# Copyright (c) 2013-2014, University Corporation for Atmospheric Research
#
# This software is distributed under a two-clause BSD license, with no
# warranties, express or implied. See the accompanying LICENSE file for
# details.
#==========================================================================

#=================================================
# Enable CTest tests.
#=================================================

enable_testing()

#=================================================
# Color output
#=================================================

option(USE_COLOR "Allow color from the build output." ON)

set(CMAKE_COLOR_MAKEFILE "${USE_COLOR}")

#=================================================
# Compiler info
#=================================================

include(Compilers)

#=================================================
# GenF90
#=================================================

option(ENABLE_GENF90
  "Use genf90.pl to regenerate out-of-date Fortran files from .in files."
  OFF)

if(ENABLE_GENF90)
  find_program(GENF90 genf90.pl)

  if(NOT GENF90)
    message(FATAL_ERROR "ENABLE_GENF90 enabled, but genf90.pl not found!")
  endif()

endif()

# Preprocessing utility functions.
include(genf90_utils)

#=================================================
# pFUnit
#=================================================

# pFUnit and its preprocessor
find_package(pFUnit REQUIRED)

# Preprocessor and driver handling.
include(pFUnit_utils)

#=================================================
# Source list and path utilities.
#=================================================

include(Sourcelist_utils)
