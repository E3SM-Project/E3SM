# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test iop remap
  CreateUnitTest(iop_remapper "iop_remapper_tests.cpp"
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test coarsening remap
  CreateUnitTest(coarsening_remapper "coarsening_remapper_tests.cpp"
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  if (EAMXX_ENABLE_EXPERIMENTAL_CODE)
    # Test refining remap (RMA version)
    CreateUnitTest(refining_remapper_rma "refining_remapper_rma_tests.cpp"
      LIBS eamxx_remap
      MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
  endif()

  # Test refining remap (P2P version)
  CreateUnitTest(refining_remapper_p2p "refining_remapper_p2p_tests.cpp"
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test vertical remap
  CreateUnitTest(vertical_remapper "vertical_remapper_tests.cpp"
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
endif()
