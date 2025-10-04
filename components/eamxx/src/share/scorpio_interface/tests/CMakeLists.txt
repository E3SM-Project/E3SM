# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  ## Test scorpio interfaces
  CreateUnitTest(scorpio_interface_tests "scorpio_interface_tests.cpp"
    LIBS eamxx_scorpio_interface LABELS "io"
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )
endif()
