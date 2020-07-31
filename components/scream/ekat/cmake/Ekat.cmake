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
                BRIEF_DOCS "Whether ekat subdir (with precision set to double) has already been processed"
                FULL_DOCS "This property is used by cmake to ensure that EKAT
                           submodule directory is only processed once (with add_subdirectory).")
define_property(GLOBAL
                PROPERTY EKAT_SINGLE_BUILT
                BRIEF_DOCS "Whether ekat subdir (with precision set to single) has already been processed"
                FULL_DOCS "This property is used by cmake to ensure that EKAT
                           submodule directory is only processed once (with add_subdirectory).")

get_property(IS_EKAT_DOUBLE_BUILT GLOBAL PROPERTY EKAT_DOUBLE_BUILT SET)
get_property(IS_EKAT_SINGLE_BUILT GLOBAL PROPERTY EKAT_SINGLE_BUILT SET)

# This macro builds EKAT with the desired precision (if not built already)
# The user can also pass an additional variable PREFIX, which will be used
# to initialize all the EKAT_BLAH variable to the value stored in ${PREFIX}_BLAH.
# If ${PREFIX}_BLAH is not set, a default will be used.
# Furthermore, we give the user the ability to override ${PREFIX}_BLAH by
# explicitly passing BLAH blahVal in the macro call.
# E.g., consider this:
#
#  buildEkat ("DOUBLE" PREFIX "MY_PROJECT" TEST_MAX_THREADS 1)
#
# It would build Ekat in double precision, set all the EKAT_XYZ var to match ${MY_PROJECT}_XYZ,
# but it would make sure that EKAT_TEST_MAX_THREADS is 1, regardless of what the
# corresponding MY_PROJECT_TEST_MAX_THREADS is
macro (buildEkat PRECISION)

  string (TOUPPER "${PRECISION}" EKAT_PRECISION)

  if (NOT "${EKAT_PRECISION}" STREQUAL "DOUBLE" AND
      NOT "${EKAT_PRECISION}" STREQUAL "SINGLE")
    message (FATAL_ERROR "Error! Unsupported precision for ekat. Valid options are: 'SINGLE', 'DOUBLE'.")
  endif ()

  # Process EKAT only if not already done
  if (NOT IS_EKAT_${EKAT_PRECISION}_BUILT)

    # First, build kokkos (if not already built), cause some of our defaults depend on it
    # Note: if kokkos is already built, this is a no-op
    include (Kokkos)

    set(options)
    set(oneValueArgs
      PREFIX
      TEST_MAX_THREADS
      TEST_THREADS_INC
      PACK_SIZE
      SMALL_PACK_SIZE
      POSSIBLY_NO_PACK
      POSSIBLY_NO_PACK_SIZE
      MPI_ERRORS_ARE_FATAL
      CONSTEXPR_ASSERT
      MIMIC_GPU
      STRICT_FP
      ENABLE_TESTS
    )
    set(multiValueArgs)
    cmake_parse_arguments(BUILD_EKAT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    # First, set EKAT_BLAH to ${${PREFIX}_BLAH}. If not set, DO set defaults
    setVars("${BUILD_EKAT_PREFIX}" TRUE)

    # Then parse the optional input, and if set, override the existing value.
    # DO NOT set defaults, or you may override something
    setVars("BUILD_EKAT" FALSE)

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
  else ()
    message ("Using Ekat with ${PRECISION} precision previously configured in this project.\n")
  endif()
endmacro()

# Set EKAT_BLAH variables from ${PREFIX}_BLAH variables. If the latter are not
# defined, use some defaults
# Note: all the EKAT variables MUST be cache variables, or EKAT will overwrite them
#       when calling set (EKAT_BLAH value TYPE CACHE "")
macro (setVars PREFIX SET_DEFAULTS)

  # Determine if this is a debug build
  string(TOLOWER "${CMAKE_BUILD_TYPE}" setVars_CMAKE_BUILD_TYPE_ci)
  if ("${setVars_CMAKE_BUILD_TYPE_ci}" STREQUAL "debug")
    set (setVars_DEBUG_BUILD TRUE)
  else ()
    set (setVars_DEBUG_BUILD FALSE)
  endif ()

  # Determine if this is a Cuda build.
  string(FIND "${KOKKOS_GMAKE_DEVICES}" "Cuda" setVars_CUDA_POS)
  if (${setVars_CUDA_POS} GREATER -1)
    set(setVars_CUDA_BUILD TRUE)
  else ()
    set(setVars_CUDA_BUILD FALSE)
  endif ()

  # Testing variables
  if (DEFINED ${PREFIX}_TEST_MAX_THREADS)
    set (EKAT_TEST_MAX_THREADS ${${PREFIX}_TEST_MAX_THREADS} CACHE STRING "")
  elseif (SET_DEFAULTS)
    set (EKAT_TEST_MAX_THREADS 1 CACHE STRING "")
  endif()

  if (DEFINED ${PREFIX}_TEST_THREAD_INC)
    set (EKAT_TEST_THREAD_INC ${${PREFIX}_TEST_THREAD_INC} CACHE STRING "")
  elseif (SET_DEFAULTS)
    set (EKAT_TEST_THREAD_INC 1 CACHE STRING "")
  endif()

  if (DEFINED ${PREFIX}_MIMIC_GPU)
    set (EKAT_MIMIC_GPU ${${PREFIX}_MIMIC_GPU} CACHE BOOL "")
  endif()

  if (DEFINED ${PREFIX}_FPMODEL)
    set (EKAT_FPMODEL ${${PREFIX}_FPMODEL} CACHE STRING "")
  elseif (setVars_DEBUG_BUILD AND SET_DEFAULTS)
    set (EKAT_FPMODEL "strict" CACHE STRING "")
  endif()

  if (DEFINED ${PREFIX}_ENABLE_TESTS)
    set (EKAT_ENABLE_TESTS ${${PREFIX}_ENABLE_TESTS} CACHE BOOL "")
  elseif (SET_DEFAULTS)
    set (EKAT_ENABLE_TESTS ON CACHE BOOL "")
  endif()

  if (DEFINED ${PREFIX}_DISABLE_TPL_WARNINGS)
    set (EKAT_DISABLE_TPL_WARNINGS ${${PREFIX}_DISABLE_TPL_WARNINGS} CACHE BOOL "")
  elseif (SET_DEFAULTS)
    set (EKAT_DISABLE_TPL_WARNINGS OFF CACHE BOOL "")
  endif()

  # Packs variables
  if (DEFINED ${PREFIX}_PACK_SIZE)
    set (EKAT_PACK_SIZE ${${PREFIX}_PACK_SIZE} CACHE STRING "")
  elseif (SET_DEFAULTS)
    if (EKAT_CUDA_BUILD)
      set (EKAT_PACK_SIZE 1 CACHE STRING "")
    else ()
      set (EKAT_PACK_SIZE 16 CACHE STRING "")
    endif()
  endif()

  if (DEFINED ${PREFIX}_SMALL_PACK_SIZE)
    set (EKAT_SMALL_PACK_SIZE ${${PREFIX}_SMALL_PACK_SIZE} CACHE STRING "")
  elseif (SET_DEFAULTS)
    set (EKAT_SMALL_PACK_SIZE ${EKAT_PACK_SIZE} CACHE STRING "")
  endif()

  if (DEFINED ${PREFIX}_PACK_CHECK_BOUNDS)
    set (EKAT_PACK_CHECK_BOUNDS ${${PREFIX}_PACK_CHECK_BOUNDS} CACHE BOOL "")
  endif()

  if (DEFINED ${PREFIX}_POSSIBLY_NO_PACK)
    set (EKAT_POSSIBLY_NO_PACK ${${PREFIX}_POSSIBLY_NO_PACK} CACHE BOOL "")
  endif()

  if (DEFINED ${PREFIX}_POSSIBLY_NO_PACK_SIZE)
    set (EKAT_POSSIBLY_NO_PACK_SIZE ${${PREFIX}_POSSIBLY_NO_PACK_SIZE} CACHE STRING "")
  elseif (SET_DEFAULTS)
    set (EKAT_POSSIBLY_NO_PACK_SIZE ${EKAT_PACK_SIZE} CACHE STRING "")
  endif()

  # MPI
  if (DEFINED ${PREFIX}_MPI_ERRORS_ARE_FATAL)
    set (EKAT_MPI_ERRORS_ARE_FATAL ${${PREFIX}_MPI_ERRORS_ARE_FATAL} CACHE BOOL "")
  endif()

  # Cleanup
  unset (setVars_CMAKE_BUILD_TYPE_ci)
  unset (setVars_DEBUG_BUILD)
  unset (setVars_CUDA_BUILD)
  unset (setVars_CUDA_POS)
endmacro ()
