include (ScreamUtils)

set (TEST_BASE_NAME p3_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LABELS p3 physics
  LIBS p3
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)  # 1h 2.5h 24h
set (RUN_T0 2021-10-12-45000)

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS p3 physics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)

# Check tendency calculation
foreach (NRANKS RANGE ${TEST_RANK_START} ${TEST_RANK_END})
  set (script ${SCREAM_BASE_DIR}/scripts/check-tendencies)
  set (fname p3_standalone_output.INSTANT.nsteps_x1.np${NRANKS}.${RUN_T0}.nc)
  set (tname p3_tend_check_np${NRANKS})
  add_test (NAME ${tname}
    COMMAND ${script} -f ${fname} -v qc T_mid -t p3_qc_tend p3_T_mid_tend
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (${tname} PROPERTIES
    LABELS "p3;physics"
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${NRANKS}_omp1)
endforeach()

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.np${TEST_RANK_END}.${RUN_T0}.nc)
  CreateBaselineTest(${TEST_BASE_NAME} ${TEST_RANK_END} ${OUT_FILE} ${FIXTURES_BASE_NAME}) 

  # Grab all bfb hash lines into a text file (only np1_omp1 test)
  CreateUnitTest (${TEST_BASE_NAME}_bfb_hash "p3_bfb_hash.cpp"
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np1_omp1
    FIXTURES_SETUP ${TEST_BASE_NAME}_gen_bfb_hash
    EXE_ARGS "--args --log-file atm.log.np1 --output-file p3_bfb_hash.txt"
    LABELS p3 baseline_gen
  )

  # Compare against bfb hashes in baseline dir
  set (SRC_FILE ${CMAKE_CURRENT_BINARY_DIR}/p3_bfb_hash.txt)
  set (TGT_FILE ${SCREAM_BASELINES_DIR}/data/p3_bfb_hash.txt)
  add_test (NAME ${TEST_BASE_NAME}_bfb_hash_baseline_cmp
    COMMAND ${CMAKE_COMMAND} -P ${SCREAM_BASE_DIR}/cmake/CompareTextFiles.cmake ${SRC_FILE} ${TGT_FILE})

  set_tests_properties (
    ${TEST_BASE_NAME}_bfb_hash_baseline_cmp PROPERTIES
    FIXTURES_REQUIRED ${TEST_BASE_NAME}_gen_bfb_hash
    LABELS "p3;baseline_cmp")

  # Store p3_bfb_hash.txt in the baseline folder
  file (APPEND ${SCREAM_TEST_OUTPUT_DIR}/baseline_list
    "${SRC_FILE}\n"
  )
endif()
