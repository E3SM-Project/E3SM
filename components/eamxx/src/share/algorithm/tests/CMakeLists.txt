# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Generate data for data interpolation tests
  CreateUnitTest(data_interpolation_setup
    SOURCES data_interpolation_setup.cpp
    LIBS eamxx_algorithm
    FIXTURES_SETUP data_interpolation_setup)

  # --- Existing time-dependent interpolation tests (split for parallelism) ---

  # Error/exception handling
  CreateUnitTest(data_interpolation_exceptions
    SOURCES data_interpolation_exceptions_tests.cpp
    LIBS eamxx_algorithm
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    FIXTURES_REQUIRED data_interpolation_setup)

  # YearlyPeriodic × Linear interpolation
  CreateUnitTest(data_interpolation_periodic_linear
    SOURCES data_interpolation_periodic_linear_tests.cpp
    LIBS eamxx_algorithm
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    FIXTURES_REQUIRED data_interpolation_setup)

  # YearlyPeriodic × Nearest interpolation
  CreateUnitTest(data_interpolation_periodic_nearest
    SOURCES data_interpolation_periodic_nearest_tests.cpp
    LIBS eamxx_algorithm
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    FIXTURES_REQUIRED data_interpolation_setup)

  # Linear timeline × Linear interpolation
  CreateUnitTest(data_interpolation_linear_hist_linear
    SOURCES data_interpolation_linear_hist_linear_tests.cpp
    LIBS eamxx_algorithm
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    FIXTURES_REQUIRED data_interpolation_setup)

  # Linear timeline × Nearest interpolation
  CreateUnitTest(data_interpolation_linear_hist_nearest
    SOURCES data_interpolation_linear_hist_nearest_tests.cpp
    LIBS eamxx_algorithm
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    FIXTURES_REQUIRED data_interpolation_setup)
  # Time interpolation
  CreateUnitTest(time_interpolation
    SOURCES time_interpolation_tests.cpp
    LIBS eamxx_algorithm
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
endif()
