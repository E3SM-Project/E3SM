# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test field manager
  CreateUnitTest(field_mgr "field_mgr_tests.cpp" LIBS eamxx_data_managers)

  ## Test single-column reader
  CreateUnitTest(scm_reader "scm_reader.cpp"
    LIBS eamxx_data_managers
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )
endif()
