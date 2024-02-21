include(ScreamUtils)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  CreateUnitTest(physics_test_data physics_test_data_unit_tests.cpp
    LIBS physics_share
    THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
endif()

if (SCREAM_ENABLE_BASELINE_TESTS)
  set(BASELINE_FILE_ARG "-b ${SCREAM_TEST_DATA_DIR}/physics_saturation.baseline")

  # The comparison test. Expects baseline to exist. All thread configurations
  # will use the same baseline.
  CreateUnitTest(physics_saturation_run_and_cmp "physics_saturation_run_and_cmp.cpp"
    LIBS physics_share
    EXE_ARGS "${BASELINE_FILE_ARG}"
    LABELS "physics")

  #
  # Use fake tests to generate shell commands to generate baselines
  #
  CreateUnitTestFromExec(physics_saturation_baseline_fake physics_saturation_run_and_cmp
    EXE_ARGS "-g ${BASELINE_FILE_ARG}"
    PROPERTIES DISABLED True)

  get_test_property(physics_saturation_baseline_fake FULL_TEST_COMMAND PHYSICS_SATURATION_GEN)

  if (PHYSICS_SATURATION_GEN STREQUAL "NOTFOUND")
    message(FATAL_ERROR "Could not get FULL_TEST_COMMAND for physics_saturation_baseline fake test")
  endif()

  separate_arguments(PHYSICS_SATURATION_GEN_ARGS UNIX_COMMAND "${PHYSICS_SATURATION_GEN}")

  add_custom_target(physics_saturation_baseline
    COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=${SCREAM_TEST_MAX_THREADS} ${PHYSICS_SATURATION_GEN_ARGS})

  add_dependencies(baseline physics_saturation_baseline)
endif()
