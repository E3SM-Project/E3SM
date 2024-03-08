include (ScreamUtils)

set (TEST_BASE_NAME spa_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LABELS spa physics
  LIBS spa
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 2 4 48)
set (ATM_TIME_STEP 864000)
set (RUN_T0 2021-10-12-45000)

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  scream/init/${EAMxx_tests_IC_FILE_72lev}
  scream/maps/map_ne4np4_to_ne2np4_mono.nc
  scream/init/spa_file_unified_and_complete_ne4_20220428.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS spa physics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX
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
