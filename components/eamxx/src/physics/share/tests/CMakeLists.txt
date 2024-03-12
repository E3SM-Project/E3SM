include(ScreamUtils)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  CreateUnitTest(physics_test_data physics_test_data_unit_tests.cpp
    LIBS physics_share
    THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
endif()

if (SCREAM_ENABLE_BASELINE_TESTS)
  if (SCREAM_BASELINES_ONLY)
    set(BASELINE_FILE_ARG "-g -b ${SCREAM_TEST_DATA_DIR}/physics_saturation.baseline")
  else()
    set(BASELINE_FILE_ARG "-b ${SCREAM_TEST_DATA_DIR}/physics_saturation.baseline")
  endif()

  # The comparison test. Expects baseline to exist. All thread configurations
  # will use the same baseline.
  CreateUnitTest(physics_saturation_run_and_cmp "physics_saturation_run_and_cmp.cpp"
    LIBS physics_share
    EXE_ARGS "${BASELINE_FILE_ARG}"
    LABELS "physics")

  # Note: the baseline_gen label is really only used if SCREAM_BASELINES_ONLY=ON, but no harm adding it
  set_tests_properties (physics_saturation_run_and_cmp PROPERTIES LABELS baseline_gen)
endif()
