# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test column ops
  CreateUnitTest(column_ops "column_ops.cpp")

  # Generate data for data interpolation test
  CreateUnitTest(data_interpolation_setup "data_interpolation_setup.cpp"
    LIBS scream_io
    FIXTURES_SETUP data_interpolation_setup)

  # Test data interpolation
  CreateUnitTest(data_interpolation "data_interpolation_tests.cpp"
    LIBS scream_io
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    FIXTURES_REQUIRED data_interpolation_setup)

  # Test common physics functions
  CreateUnitTest(common_physics "common_physics_functions_tests.cpp")

  # Test atmosphere processes
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/atm_process_tests_named_procs.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/atm_process_tests_named_procs.yaml COPYONLY)

  CreateUnitTest(atm_proc "atm_process_tests.cpp")

endif()
