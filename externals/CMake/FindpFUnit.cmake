# Find module for pFUnit
#
# For this module to work, either the pFUnit parser must be discoverable
# (e.g. in the user's PATH), or else the environment variable "PFUNIT" must
# be defined, and point to the root directory for the PFUNIT installation.
#
# This module sets some typical variables:
# PFUNIT_FOUND
# PFUNIT_LIBRARY(/LIBRARIES)
# PFUNIT_INCLUDE_DIR(/DIRS)
#
# The module also sets:
# PFUNIT_DRIVER - Path to the pFUnit driver source.
# PFUNIT_MODULE_DIR - Directory containing pFUnit's module files.
# PFUNIT_PARSER - Path to pFUnitParser.py (the preprocessor).

#==========================================================================
# Copyright (c) 2013-2014, University Corporation for Atmospheric Research
#
# This software is distributed under a two-clause BSD license, with no
# warranties, express or implied. See the accompanying LICENSE file for
# details.
#==========================================================================

include(FindPackageHandleStandardArgs)

find_program(PFUNIT_PARSER pFUnitParser.py
  HINTS $ENV{PFUNIT}/bin)

string(REGEX REPLACE "bin/pFUnitParser\\.py\$" ""
  pfunit_directory ${PFUNIT_PARSER})

find_library(PFUNIT_LIBRARY pfunit
  HINTS ${pfunit_directory}/lib)

find_path(PFUNIT_INCLUDE_DIR driver.F90
  HINTS ${pfunit_directory}/include)

set(PFUNIT_DRIVER ${PFUNIT_INCLUDE_DIR}/driver.F90)

find_path(PFUNIT_MODULE_DIR NAMES pfunit.mod PFUNIT.MOD
  HINTS ${pfunit_directory}/include ${pfunit_directory}/mod)

set(PFUNIT_LIBRARIES ${PFUNIT_LIBRARY})
set(PFUNIT_INCLUDE_DIRS ${PFUNIT_INCLUDE_DIR} ${PFUNIT_MODULE_DIR})

# Handle QUIETLY and REQUIRED.
find_package_handle_standard_args(pFUnit DEFAULT_MSG
  PFUNIT_LIBRARY PFUNIT_INCLUDE_DIR PFUNIT_MODULE_DIR PFUNIT_PARSER)

mark_as_advanced(PFUNIT_INCLUDE_DIR PFUNIT_LIBRARY PFUNIT_MODULE_DIR
  PFUNIT_PARSER)
