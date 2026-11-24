# Homme tests. See README.md in this folder for an overview.
#
#   execs/  : the executables used by the tests (each one is built in <build-dir>/test/execs/<name>)
#   unit/   : unit tests of the C++ code (often comparing against the F90 implementation)
#   system/ : tests that run a Homme executable with some namelists, and check the output
#             (against baselines, against the F90/C++ counterpart, ...)
#   cases/  : setups for hand-run simulations (e.g., DCMIP cases). They are not ctest tests.
#   data/   : input files used by tests and cases (vertical coordinates, meshes)
#   cmake/  : the CMake utilities used to create the tests, and the scripts that ctest runs
#
# Targets:
#   test-execs: build all the executables needed by the tests (and cprnc)
#   check     : run all the tests (same as running ctest)
#   baseline  : run the tests that generate baselines (ctest -L baseline_gen), and copy
#               their outputs in <build-dir>/tests/baseline
# NOTE: CIME's HOMME test uses the test-execs, check, and baseline targets.

set (HOMME_TEST_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)
set (HOMME_SYSTEM_TESTS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/system)

include (HommeTesting)
homme_setup_testing()

add_custom_target (test-execs)
if (NOT BUILD_HOMME_WITHOUT_PIOLIBRARY)
  # The tests need cprnc (a no-op if cprnc is a prebuilt binary, found via CPRNC_DIR)
  add_dependencies (test-execs cprnc)
endif()

add_custom_target (check
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_dependencies (check test-execs)

add_custom_target (baseline
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L baseline_gen
  COMMAND ${CMAKE_COMMAND} -P ${HOMME_COPY_BASELINES_SCRIPT}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_dependencies (baseline test-execs)

add_subdirectory (execs)
add_subdirectory (unit)
add_subdirectory (system)
add_subdirectory (cases)

homme_write_copy_baselines_script()
