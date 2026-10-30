# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test iop remap
  CreateUnitTest(iop_remapper
    SOURCES iop_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  # Test coarsening remap
  CreateUnitTest(coarsening_remapper
    SOURCES coarsening_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  # Test refining remap
  CreateUnitTest(refining_remapper
    SOURCES refining_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  # Test remap from a grid to its unique grid
  CreateUnitTest(unique_grid_remapper
    SOURCES unique_grid_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  # Test vertical remap
  CreateUnitTest(vertical_remapper
    SOURCES vertical_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )
endif()
