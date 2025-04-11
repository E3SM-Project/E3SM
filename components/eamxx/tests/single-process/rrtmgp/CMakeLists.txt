INCLUDE (ScreamUtils)

set (TEST_BASE_NAME rrtmgp_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})

set(YAKL_LIB_NAME "")
if (SCREAM_RRTMGP_ENABLE_YAKL)
  set(YAKL_LIB_NAME "yakl")
endif()

if (SCREAM_ENABLE_BASELINE_TESTS AND NOT SCREAM_ONLY_GENERATE_BASELINES)
  # Unit test to compare against raw rrtmgp output
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_unit.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_unit.yaml)
  CreateUnitTest(${TEST_BASE_NAME}_unit rrtmgp_standalone_unit.cpp
    LABELS rrtmgp physics driver
    LIBS scream_rrtmgp rrtmgp scream_control ${YAKL_LIB_NAME} diagnostics rrtmgp_test_utils
    EXE_ARGS "--args --rrtmgp_inputfile ${SCREAM_DATA_DIR}/init/rrtmgp-allsky.nc --rrtmgp_baseline ${SCREAM_BASELINES_DIR}/data/rrtmgp-allsky-baseline.nc"
  )
endif()

## Create rrtmgp stand alone executable
CreateUnitTestExec(${TEST_BASE_NAME} "rrtmgp_standalone.cpp"
  LIBS scream_rrtmgp rrtmgp scream_control ${YAKL_LIB_NAME} diagnostics
)

# The RRTMGP stand-alone test that runs multi-step
# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Test non-chunked version (sweep multiple ranks)
set (SUFFIX "_not_chunked")
set (COL_CHUNK_SIZE 1000)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                ${CMAKE_CURRENT_BINARY_DIR}/output_not_chunked.yaml)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                ${CMAKE_CURRENT_BINARY_DIR}/input_not_chunked.yaml)
CreateUnitTestFromExec(
  ${TEST_BASE_NAME}_not_chunked ${TEST_BASE_NAME}
  LABELS rrtmgp physics driver
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  EXE_ARGS "--args -inputfile input_not_chunked.yaml"
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}_not_chunked
)

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}_not_chunked
  FILE_META_NAME ${TEST_BASE_NAME}_output_not_chunked.INSTANT.nsteps_x${NUM_STEPS}.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS rrtmgp physics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_not_chunked_npMPIRANKS_omp1
)

## Test chunked version (only for ${TEST_RANK_END}) and compare against non-chunked
set (SUFFIX "_chunked")
math (EXPR COL_PER_RANK "218 / ${TEST_RANK_END}")
math (EXPR COL_CHUNK_SIZE "${COL_PER_RANK} / 2")
if (COL_CHUNK_SIZE LESS 1)
  message (FATAL_ERROR "Error! Chunk size for rrtmgp unit test is less than 1.")
endif()

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                ${CMAKE_CURRENT_BINARY_DIR}/input_chunked.yaml)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                ${CMAKE_CURRENT_BINARY_DIR}/output_chunked.yaml)
CreateUnitTestFromExec(
    ${TEST_BASE_NAME}_chunked ${TEST_BASE_NAME}
    LABELS rrtmgp physics driver
    MPI_RANKS ${TEST_RANK_END}
    EXE_ARGS "--args --inputfile=input_chunked.yaml"
    FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}_chunked
    PROPERTIES PASS_REGULAR_EXPRESSION "(beg.end: 0, ${COL_CHUNK_SIZE})"
)

CompareNCFiles(
  TEST_NAME ${TEST_BASE_NAME}_chunked_vs_not_chunked
  SRC_FILE ${TEST_BASE_NAME}_output_chunked.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.${RUN_T0}.nc
  TGT_FILE ${TEST_BASE_NAME}_output_not_chunked.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.${RUN_T0}.nc
  LABELS rrtmgp physics
  FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_chunked_np${TEST_RANK_END}_omp1
                    ${FIXTURES_BASE_NAME}_not_chunked_np${TEST_RANK_END}_omp1)

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX,
  #       and that chunked is bfb with not_chunked
  set (OUT_FILE ${TEST_BASE_NAME}_output_chunked.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.${RUN_T0}.nc)
  CreateBaselineTest(${TEST_BASE_NAME}_chunked ${TEST_RANK_END} ${OUT_FILE} ${FIXTURES_BASE_NAME}_chunked)
endif()
