include (ScreamUtils)

# This function compares an output file created by an AD-driven test with its
# corresponding one in the baselines dir. The inputs (all REQUIRED) are
#
# - TEST_BASE_NAME: the base name of the test that generated the file; will be used also
#                   as the prefix for the baseline_cmp test name
# - GEN_TEST_NRANKS: the number of ranks used in the test that generated the file
# - OUT_FILE: the name of the output nc file
# - FIXTURES_BASE_NAME: base name of the FIXTURES_SETUP property of the test that generated the file.
#                       The string _np${GEN_TEST_NRANKS}_omp1 will be attached to this

function (CreateBaselineTest TEST_BASE_NAME GEN_TEST_NRANKS OUT_FILE FIXTURES_BASE_NAME)
  # Get names of src and tgt files
  set (SRC_FILE ${CMAKE_CURRENT_BINARY_DIR}/${OUT_FILE})
  set (TGT_FILE ${SCREAM_BASELINES_DIR}/data/${OUT_FILE})

  # Add comparison test using the CprncTest.cmake scritp shipped by EAMxx
  add_test (
    NAME ${TEST_BASE_NAME}_baseline_cmp
    COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(
    ${TEST_BASE_NAME}_baseline_cmp
    PROPERTIES
    LABELS baseline_cmp
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${TEST_RANK_END}_omp1)

  # Add the test that generated the baseline to the baseline_gen label, so scripts/test-all-eamxx
  # can run it when it has to generate baselines
  if (TEST ${TEST_BASE_NAME}_np${GEN_TEST_NRANKS})
    set (GEN_TEST_FULL_NAME ${TEST_BASE_NAME}_np${GEN_TEST_NRANKS})
  elseif(TEST ${TEST_BASE_NAME})
    set (GEN_TEST_FULL_NAME ${TEST_BASE_NAME})
  else()
    string (CONCAT msg
      "Could not find the test that generated the output file\n"
      "  TEST_BASE_NAME: ${TEST_BASE_NAME}\n"
      "  OUT_FILE      : ${OUT_FILE}\n")
    message ("${msg}")
    message (FATAL_ERROR "Aborting...")
  endif()
  set_tests_properties (${GEN_TEST_FULL_NAME} PROPERTIES LABELS baseline_gen)

  # test-all-eamxx will read this file to get the list of baseline nc files to
  # copy into the baseline dir
  file (APPEND ${SCREAM_TEST_OUTPUT_DIR}/baseline_list
    "${SRC_FILE}\n"
  )

endfunction()

# Some tests for checking that certain testing infrastructures work
add_subdirectory(generic)

if (NOT DEFINED ENV{SCREAM_FAKE_ONLY})
  # memcheck builds (and also coverage ones) can just run the max ranks, since they
  # do only need to perform checks on the code itself, rather than the model it represents
  SetVarDependingOnTestSize(TEST_RANK_START ${SCREAM_TEST_MAX_RANKS} 1 1)
  set(TEST_RANK_END ${SCREAM_TEST_MAX_RANKS})

  # Initial condition files used in the tests
  set(EAMxx_tests_IC_FILE_72lev         "screami_unit_tests_ne2np4L72_20220822.nc")
  set(EAMxx_tests_IC_FILE_128lev        "screami_unit_tests_ne2np4L128_20220822.nc")
  set(EAMxx_tests_TOPO_FILE             "USGS-gtopo30_ne2np4pg2_x6t_20230331.nc")
  set(EAMxx_tests_IC_FILE_MAM4xx_72lev  "screami_unit_tests_mam4xx_ne2np4L72_20240329.nc")

  # Testing individual atm processes
  add_subdirectory(single-process)

  # Testing multiple atm processes coupled together.
  # Some compute-sanitizer tests time out with these
  # larger multiprocess tests, so disable in that case.
  if (NOT EKAT_ENABLE_COMPUTE_SANITIZER)
    add_subdirectory(multi-process)
  endif()
  
  if (EAMXX_ENABLE_PYSCREAM)
    add_subdirectory(python)
  endif()
endif()
