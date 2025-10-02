# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test point grid
  CreateUnitTest(point_grid "point_grid_tests.cpp"
    LIBS eamxx_grid
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test se grid
  CreateUnitTest(se_grid "se_grid_tests.cpp"
    LIBS eamxx_grid
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test grid gid-related methods
  CreateUnitTest(grid "grid_tests.cpp"
    LIBS eamxx_grid
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test grid import-export
  CreateUnitTest(grid_imp_exp "grid_import_export_tests.cpp"
    LIBS eamxx_grid
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
endif()
