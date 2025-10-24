# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Generate data for data interpolation test
  CreateUnitTest(data_interpolation_setup "data_interpolation_setup.cpp"
    LIBS eamxx_algorithm
    FIXTURES_SETUP data_interpolation_setup)

  # Test data interpolation
  CreateUnitTest(data_interpolation "data_interpolation_tests.cpp"
    LIBS eamxx_algorithm
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    FIXTURES_REQUIRED data_interpolation_setup)

  # Time interpolation
  CreateUnitTest(time_interpolation "time_interpolation_tests.cpp"
    LIBS eamxx_algorithm
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
endif()
