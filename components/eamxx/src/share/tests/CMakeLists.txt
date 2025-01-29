# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test utils
  CreateUnitTest(utils "utils_tests.cpp")

  # Test column ops
  CreateUnitTest(column_ops "column_ops.cpp")

  # Test fields
  CreateUnitTest(field "field_tests.cpp")

  # Test fields
  CreateUnitTest(subfield "subfield_tests.cpp")

  # Test field utils
  CreateUnitTest(field_utils "field_utils.cpp"
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test property checks
  CreateUnitTest(property_checks "property_checks.cpp")

  # Test grids
  CreateUnitTest(grid "grid_tests.cpp"
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test grid import-export
  CreateUnitTest(grid_imp_exp "grid_import_export_tests.cpp"
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test iop remap
  CreateUnitTest(iop_remapper "iop_remapper_tests.cpp"
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test coarsening remap
  CreateUnitTest(coarsening_remapper "coarsening_remapper_tests.cpp"
    LIBS scream_io
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  if (EAMXX_ENABLE_EXPERIMENTAL_CODE)
    # Test refining remap (RMA version)
    CreateUnitTest(refining_remapper_rma "refining_remapper_rma_tests.cpp"
      LIBS scream_io
      MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})
  endif()

  # Test refining remap (P2P version)
  CreateUnitTest(refining_remapper_p2p "refining_remapper_p2p_tests.cpp"
    LIBS scream_io
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  # Test vertical remap
  CreateUnitTest(vertical_remapper "vertical_remapper_tests.cpp"
    LIBS scream_io
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

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
