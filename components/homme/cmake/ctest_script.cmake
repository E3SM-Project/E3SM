cmake_minimum_required(VERSION 3.16)

# ---------------------------------------------------------------------------
# Input variables (all may be overridden via -D on the ctest command line)
# ---------------------------------------------------------------------------
if (NOT DEFINED HOMME_ROOT)
  get_filename_component(HOMME_ROOT "${CMAKE_CURRENT_LIST_DIR}" DIRECTORY)
endif()

if (NOT DEFINED BUILD_WORK_DIR)
  set(BUILD_WORK_DIR "${HOMME_ROOT}/ctest-build")
endif()

if (NOT DEFINED MACHINE_FILE)
  message(FATAL_ERROR "MACHINE_FILE is required (e.g., -DMACHINE_FILE=pm-cpu-bfb)")
endif()

if (NOT DEFINED BFB_TESTING)
  set(BFB_TESTING FALSE)
endif()

if (NOT DEFINED BUILD_PARALLEL_LEVEL)
  set(BUILD_PARALLEL_LEVEL 8)
endif()

if (NOT DEFINED RUN_PARALLEL_LEVEL)
  set(RUN_PARALLEL_LEVEL 4)
endif()

if (NOT DEFINED DEBUG)
  set(DEBUG FALSE)
endif()

if (NOT DEFINED SUBMIT)
  set(SUBMIT FALSE)
endif()

if (NOT DEFINED dashboard_model)
  set(dashboard_model Experimental)
endif()

if (NOT DEFINED dashboard_track)
  set(dashboard_track HOMME)
endif()

# ---------------------------------------------------------------------------
# Derived settings
# ---------------------------------------------------------------------------
# MACHINE_FILE is the name (without .cmake) of a file in cmake/machineFiles
set(_machine_file "${HOMME_ROOT}/cmake/machineFiles/${MACHINE_FILE}.cmake")
if (NOT EXISTS "${_machine_file}")
  message(FATAL_ERROR "Machine file not found: ${_machine_file}")
endif()

if (NO_BASELINE)
  # Tests will not be compared against baselines (only against each other, e.g., f90 vs cxx)
  if (GENERATE)
    message (FATAL_ERROR "Cannot generate baselines when NO_BASELINE is TRUE")
  endif()
elseif (NOT DEFINED BASELINE_DIR)
  message (FATAL_ERROR "BASELINE_DIR was not specified")
endif()

# ---------------------------------------------------------------------------
# CDash / dashboard setup
# ---------------------------------------------------------------------------
set(CTEST_SOURCE_DIRECTORY "${HOMME_ROOT}")
set(CTEST_BINARY_DIRECTORY "${BUILD_WORK_DIR}")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
if (BFB_TESTING)
  set(CTEST_BUILD_NAME "homme")
else()
  set(CTEST_BUILD_NAME "homme-bfb")
endif()

if (DEFINED ENV{HOSTNAME})
  set(CTEST_SITE "$ENV{HOSTNAME}")
else()
  cmake_host_system_information(RESULT CTEST_SITE QUERY HOSTNAME)
endif()

file(MAKE_DIRECTORY "${BUILD_WORK_DIR}")

ctest_start(${dashboard_model} TRACK ${dashboard_track})

# ---------------------------------------------------------------------------
# Configure
# ---------------------------------------------------------------------------
set(_configure_opts "-C${_machine_file}")
if (NO_BASELINE)
  list(APPEND _configure_opts "-DHOMME_SKIP_BASELINE_COMPARE=ON")
else()
  list(APPEND _configure_opts "-DHOMME_BASELINE_DIR=${BASELINE_DIR}")
endif()
list(APPEND _configure_opts "-DCMAKE_VERBOSE_MAKEFILE=ON")

if (DEBUG)
  list(APPEND _configure_opts "-DCMAKE_BUILD_TYPE=Debug")
endif()
if (CPRNC_DIR)
  list(APPEND _configure_opts "-DCPRNC_DIR=${CPRNC_DIR}")
endif()
foreach (_cm_opt IN LISTS EXTRA_CMAKE_OPTIONS)
  list(APPEND _configure_opts "${_cm_opt}")
endforeach()

ctest_configure(OPTIONS "${_configure_opts}" RETURN_VALUE _cfg_rv)

set(_had_failure FALSE)
if (_cfg_rv)
  set(_had_failure TRUE)
else()
  # -------------------------------------------------------------------------
  # Build
  # -------------------------------------------------------------------------
  if (GENERATE)
    # For generate, Homme uses 'make baseline', rather than rely on ctest labels.
    # The target depends on all the test executables, so it builds them first, then runs the tests.
    message(STATUS "HOMME ctest script: generate baselines")
    set(_build_target "baseline")
  else()
    set(_build_target "test-execs")
  endif()

  # With parallel make, --output-sync=target keeps each target's output together. Otherwise, interleaved
  # compiler output can confuse ctest's error/warning parsing, producing bogus "errors" in the build log.
  ctest_build(
    TARGET "${_build_target}"
    FLAGS "-j${BUILD_PARALLEL_LEVEL} --output-sync=target"
    RETURN_VALUE _bld_rv
  )
  if (_bld_rv)
    set(_had_failure TRUE)
  elseif (NOT GENERATE)
    message(STATUS "HOMME ctest script: compare against baselines")
    ctest_test(
      PARALLEL_LEVEL ${RUN_PARALLEL_LEVEL}
      RETURN_VALUE _test_rv
    )
    if (_test_rv)
      set(_had_failure TRUE)
    endif()
  endif()
endif()

# ---------------------------------------------------------------------------
# Upload logs and submit to CDash
# ---------------------------------------------------------------------------
if (SUBMIT)
  ctest_submit(RETRY_COUNT 10 RETRY_DELAY 60)
endif()

# ---------------------------------------------------------------------------
# Ensure the overall job returns an error code to the calling process
# ---------------------------------------------------------------------------
if (_had_failure)
  message(FATAL_ERROR "One or more HOMME ctest steps failed")
endif()
