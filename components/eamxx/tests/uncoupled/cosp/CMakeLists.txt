include (ScreamUtils)

set (TEST_BASE_NAME cosp_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LABELS cosp physics
  LIBS eamxx_cosp
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/USGS-gtopo30_ne4np4pg2_16x_converted.c20200527.nc)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: for other tests we do np1-vs-npX bfb tests, which is why one is enough.
  #       COSP, is not bfb w.r.t. num ranks though. Still, we only check 1 output files against baselines
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.np${TEST_RANK_END}.${RUN_T0}.nc)
  set (SRC_FILE ${CMAKE_CURRENT_BINARY_DIR}/${OUT_FILE})
  set (TGT_FILE ${SCREAM_BASELINES_DIR}/data/${OUT_FILE})
  add_test (
    NAME ${TEST_BASE_NAME}_baseline_cmp
    COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(
    ${TEST_BASE_NAME}_baseline_cmp
    PROPERTIES
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${TEST_RANK_END}_omp1)

  # The np${TEST_RANK_END} test is used for baselines tests, so add it to the baseline label
  if (TEST_RANK_END GREATER 1)
    # ECUT only adds _npX if we have more than one value of X, or if X>1
    set (TEST_SUFFIX _np${TEST_RANK_END})
  endif()
  set_tests_properties (${TEST_BASE_NAME}${TEST_SUFFIX} PROPERTIES LABELS baseline)

  # test-all-scream will read this file to get the list of baseline nc files to
  # copy into the baseline dir
  file (APPEND ${SCREAM_TEST_OUTPUT_DIR}/baseline_list
    "${SRC_FILE}\n"
  )
endif()
