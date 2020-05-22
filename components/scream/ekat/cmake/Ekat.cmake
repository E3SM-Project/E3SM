# This cmake utility allows you to build ekat just by adding
#
#   Include(Ekat)
#
# to your CMakeLists.txt. The advantage over manually adding the
# ekat subdirectory is that this script can check whether ekat was
# already built with the desired configuration, and avoid rebuilding
# it. For instance, say some of your subprojects need ekat in single
# precision, and some need it in double precision. We would need
# something like this
#
#  set (EKAT_DOUBLE_PRECISION ON)
#  add_subdirectory(path/to/ekat some/binary/folder/ekat_dp)
#
#  set (EKAT_DOUBLE_PRECISION OFF)
#  add_subdirectory(path/to/ekat some/binary/folder/ekat_dp)
#
# However, before adding the subdir, you need to make sure that
# the binary folder is 'available', that is, nobody has yet tried to
# build ekat with that precision. This script can do that for you

# Define global properties for EKAT_DOUBLE_BUILT and EKAT_SINGLE_BUILT
# to ensure ekat (with the desired precision) is built only once
define_property(GLOBAL
                PROPERTY EKAT_DOUBLE_BUILT
                BRIEF_DOCS "Wheter ekat subdir (with precision set to double) has already been processed"
                FULL_DOCS "This property is used by cmake to ensure that EKAT
                           submodule directory is only processed once (with add_subdirectory).")
define_property(GLOBAL
                PROPERTY EKAT_SINGLE_BUILT
                BRIEF_DOCS "Wheter ekat subdir (with precision set to single) has already been processed"
                FULL_DOCS "This property is used by cmake to ensure that EKAT
                           submodule directory is only processed once (with add_subdirectory).")

get_property(IS_EKAT_DOUBLE_BUILT GLOBAL PROPERTY EKAT_DOUBLE_BUILT SET)
get_property(IS_EKAT_SINGLE_BUILT GLOBAL PROPERTY EKAT_SINGLE_BUILT SET)

# Set precision to double if not requested
if (NOT EKAT_PRECISION)
  set (EKAT_PRECISION "DOUBLE")
endif()

if (NOT "${EKAT_PRECISION}" STREQUAL "DOUBLE" AND
    NOT "${EKAT_PRECISION}" STREQUAL "SINGLE")
  message (FATAL_ERROR "Error! Unsupported precision for ekat. Valid options are: 'SINGLE', 'DOUBLE'.")
endif ()

# Process EKAT, if not already done
if (NOT IS_EKAT_${EKAT_PRECISION}_BUILT)

  # Note: we need to set EKAT_DOUBLE_PRECISION in the cache, or ekat will overwrite it
  #       when calling  option(EKAT_DOUBLE_PRECISION ...)
  if ("${EKAT_PRECISION}" STREQUAL "DOUBLE")
    set (EKAT_DOUBLE_PRECISION ON  CACHE BOOL "Whether EKAT will be built in double precision")
  else ()
    set (EKAT_DOUBLE_PRECISION OFF CACHE BOOL "Whether EKAT will be built in double precision")
  endif()

  if (NOT EKAT_SOURCE_DIR)
    message (FATAL_ERROR "Error! Please, specify path to EKAT in EKAT_SOURCE_DIR.\n")
  endif()
  if (NOT EXISTS ${EKAT_SOURCE_DIR}/ekat_config.h.in)
    message (FATAL_ERROR "Error! Something is wrong with EKAT_SOURCE_DIR."
                         "       No 'ekat_config.h.in' file was found in '${EKAT_SOURCE_DIR}'")
  endif()
  add_subdirectory (${EKAT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/externals/ekat_${EKAT_PRECISION})

  # Make sure that future includes of this script don't rebuild the same configuration
  set_property(GLOBAL PROPERTY EKAT_${EKAT_PRECISION}_BUILT TRUE)
endif()
