include(ScreamUtils)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  CreateUnitTest(physics_test_data
    SOURCES physics_test_data_unit_tests.cpp
    LIBS eamxx_physics_share
    THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
  )

  # Test common physics functions
  CreateUnitTest(common_physics
    SOURCES common_physics_functions_tests.cpp
    LIBS eamxx_physics_share
  )
endif()

CreateUnitTest(physics_saturation_run_and_cmp
  SOURCES physics_saturation_run_and_cmp.cpp
  LIBS eamxx_physics_share
  EXE_ARGS "${SCREAM_BASELINE_FILE_ARG}"
  THREADS ${SCREAM_BASELINE_TEST_THREADS}
  LABELS saturation physics baseline_gen baseline_cmp
  )
