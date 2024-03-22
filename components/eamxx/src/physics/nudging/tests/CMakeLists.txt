if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  CreateUnitTest (nudging_data "create_nudging_data.cpp"
    LIBS scream_io
    FIXTURES_SETUP nudging_data
  )
  CreateUnitTest (nudging_map_files "create_map_file.cpp"
    LIBS scream_io
    FIXTURES_SETUP nudging_map_files
  )
  CreateUnitTest(nudging_tests "nudging_tests.cpp"
    LIBS nudging scream_io
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    LABELS physics nudging
    FIXTURES_REQUIRED "nudging_map_files;nudging_data"
  )

endif()
