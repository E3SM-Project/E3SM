# Homme tests: build the test executables, and create the ctest tests.
#
# The tests are created with the functions in cmake/HommeTesting.cmake. Each system test X
# runs in ${CMAKE_BINARY_DIR}/tests/X. Use 'ctest -N -L <label>' to list tests with a given label
# (e.g., system, unit, baseline_gen, baseline_cmp, cxx_vs_f90, theta-l_kokkos, ...).
#
# Targets:
#   test-execs: build all the executables needed by the tests (and cprnc)
#   check     : run all the tests (same as running ctest)
#   baseline  : run the tests that generate baselines (ctest -L baseline_gen), and copy
#               their outputs in ${CMAKE_BINARY_DIR}/tests/baseline
# NOTE: CIME's HOMME test uses the test-execs, check, and baseline targets.

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

# Add the test exec subdirs for the sweqx executable
IF(${BUILD_HOMME_SWEQX})
  ADD_SUBDIRECTORY(swtcA)
  ADD_SUBDIRECTORY(swtcB)
  ADD_SUBDIRECTORY(swtcC)
ENDIF()


# Add the test exec subdirs for the preqx executable
IF(${BUILD_HOMME_PREQX})
  ADD_SUBDIRECTORY(baroC)
  ADD_SUBDIRECTORY(baroCam)
  ADD_SUBDIRECTORY(preqx-nlev10-native)
  ADD_SUBDIRECTORY(preqx-nlev20-interp)
  ADD_SUBDIRECTORY(preqx-nlev20-native)
  ADD_SUBDIRECTORY(preqx-nlev30-interp)
  ADD_SUBDIRECTORY(preqx-nlev30-native)
  ADD_SUBDIRECTORY(preqx-nlev40-interp)
  ADD_SUBDIRECTORY(preqx-nlev60-interp)
  ADD_SUBDIRECTORY(preqx-nlev60-native)
#  ADD_SUBDIRECTORY(preqx-nlev100-native)
#  ADD_SUBDIRECTORY(preqx-nlev140-native)
#  ADD_SUBDIRECTORY(preqx-nlev200-native)
ENDIF()

# Add the test exec subdirs for the prim executable
IF(BUILD_HOMME_TOOL)
  ADD_SUBDIRECTORY(tool-nlev26)
ENDIF()

IF(${BUILD_HOMME_THETA})
  ADD_SUBDIRECTORY(theta-l-nlev20)
  ADD_SUBDIRECTORY(theta-l-nlev30)
  ADD_SUBDIRECTORY(theta-l-nlev40)
  ADD_SUBDIRECTORY(theta-l-nlev60)
  ADD_SUBDIRECTORY(theta-l-nlev10-native)
  ADD_SUBDIRECTORY(theta-l-nlev20-native)
#  ADD_SUBDIRECTORY(theta-l-nlev60-native)
#  ADD_SUBDIRECTORY(theta-l-nlev30-native)
#  ADD_SUBDIRECTORY(theta-l-nlev64-native)
#  ADD_SUBDIRECTORY(theta-l-nlev100-native)
  ADD_SUBDIRECTORY(theta-l-nlev128-native)
#  ADD_SUBDIRECTORY(theta-l-nlev140-native)
#  ADD_SUBDIRECTORY(theta-l-nlev150-native)
#  ADD_SUBDIRECTORY(theta-l-nlev200-native)
#  ADD_SUBDIRECTORY(theta-l-nlev256-native)
#  ADD_SUBDIRECTORY(theta-l-nlev300-native)
  # Special test for -DHOMME_WITHOUT_PIOLIBRARY.
  ADD_SUBDIRECTORY(stt)
ENDIF()

IF (${BUILD_HOMME_THETA_KOKKOS})
  ADD_SUBDIRECTORY(theta-l-nlev30-kokkos)
  ADD_SUBDIRECTORY(theta-l-nlev72-kokkos)

  ADD_SUBDIRECTORY(theta-nlev128-kokkos)

  # For doubly periodic tests.
  ADD_SUBDIRECTORY(theta-l-nlev20-native-kokkos)
  ADD_SUBDIRECTORY(theta-l-nlev10-native-kokkos)
  ADD_SUBDIRECTORY(theta-l-nlev128-native-kokkos)

  # This one we build it only if both theta and theta_kokkos are requested,
  # since it's just needed for perf/correctness comparisons
  IF (BUILD_HOMME_THETA)
    ADD_SUBDIRECTORY(theta-nlev128)
  ENDIF()

  ADD_SUBDIRECTORY(thetal_kokkos_ut)

ENDIF()

# Test some stuff in the shared folder
ADD_SUBDIRECTORY (share_kokkos_ut)

IF (${BUILD_HOMME_PREQX_KOKKOS})
  ADD_SUBDIRECTORY(preqx-nlev26-kokkos)
  ADD_SUBDIRECTORY(preqx-nlev72-kokkos)
  IF (HOMMEXX_BFB_TESTING)
    ADD_SUBDIRECTORY(preqx-nlev26)
    ADD_SUBDIRECTORY(preqx-nlev72)
    ADD_SUBDIRECTORY(preqx_kokkos_ut)
  ENDIF ()
    # For doubly periodic tests.
    ADD_SUBDIRECTORY(preqx-nlev10-native-kokkos)
ENDIF ()

# Read the test-list.cmake file to get the HOMME_TESTS list
include (${HOMME_SOURCE_DIR}/test/reg_test/run_tests/test-list.cmake)
createTests(HOMME_TESTS)

if (HOMMEXX_BFB_TESTING AND BUILD_HOMME_PREQX_KOKKOS)
  message (STATUS "Generating preqx tests no more intensive than profile ${HOMME_TESTING_PROFILE}.")
  make_profiles_up_to(${HOMME_TESTING_PROFILE} profiles)
  foreach (p ${profiles})
    message (STATUS "For profile ${p}:")
    createTestsWithProfile(HOMME_PREQX_TESTS_WITH_PROFILE ${p})
    CREATE_CXX_VS_F90_TESTS_WITH_PROFILE(PREQX_COMPARE_F_C_TEST ${p})
  endforeach ()
endif()
if (HOMMEXX_BFB_TESTING AND BUILD_HOMME_THETA_KOKKOS)
  message (STATUS "Generating theta tests no more intensive than profile ${HOMME_TESTING_PROFILE}.")
  make_profiles_up_to(${HOMME_TESTING_PROFILE} profiles)
  foreach (p ${profiles})
    message (STATUS "For profile ${p}:")
    createTestsWithProfile(HOMME_THETA_TESTS_WITH_PROFILE ${p})
    CREATE_CXX_VS_F90_TESTS_WITH_PROFILE(THETA_COMPARE_F_C_TEST ${p})
  endforeach ()
endif()

homme_write_copy_baselines_script()
