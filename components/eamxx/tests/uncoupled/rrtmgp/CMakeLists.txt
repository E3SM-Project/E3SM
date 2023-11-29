INCLUDE (ScreamUtils)

# Test atmosphere processes
if (NOT SCREAM_BASELINES_ONLY)
  # Ensure test input files are present in the data dir
  GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})

  set (TEST_BASE_NAME rrtmgp_standalone)
  set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

  # Unit test to compare against raw rrtmgp output
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_unit.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_unit.yaml)
  CreateUnitTest(${TEST_BASE_NAME}_unit rrtmgp_standalone_unit.cpp
    LABELS rrtmgp physics driver
    LIBS scream_rrtmgp rrtmgp scream_control yakl diagnostics rrtmgp_test_utils
    EXE_ARGS "--ekat-test-params rrtmgp_inputfile=${SCREAM_DATA_DIR}/init/rrtmgp-allsky.nc,rrtmgp_baseline=${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc"
  )
  # This test needs the allsky baselines file
  add_dependencies (rrtmgp_standalone_unit rrtmgp_allsky_baseline.nc)

  ## Create rrtmgp stand alone executable
  CreateUnitTestExec(${TEST_BASE_NAME} "rrtmgp_standalone.cpp"
    LIBS scream_rrtmgp rrtmgp scream_control yakl diagnostics
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
    EXE_ARGS "--ekat-test-params inputfile=input_not_chunked.yaml"
    FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}_not_chunked
  )

  # Compare output files produced by npX tests, to ensure they are bfb
  include (CompareNCFilesFamily)

  CreateRange(MpiRanks ${TEST_RANK_START} ${TEST_RANK_END})

  CompareNCFilesFamily (
    TEST_META_NAME ${TEST_BASE_NAME}_not_chunked_npMPIRANKS_vs_np1
    FILE_META_NAME ${TEST_BASE_NAME}_output_not_chunked.INSTANT.nsteps_x${NUM_STEPS}.npMPIRANKS.${RUN_T0}.nc
    MAGIC_STRING "MPIRANKS"
    MAGIC_VALUES "${MpiRanks}"
    LABELS rrtmgp physics PEM
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_not_chunked_np1_omp1
                      ${FIXTURES_BASE_NAME}_not_chunked_npMPIRANKS_omp1)


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
      EXE_ARGS "--ekat-test-params inputfile=input_chunked.yaml"
      FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}_chunked
      PROPERTIES PASS_REGULAR_EXPRESSION "(beg.end: 0, ${COL_CHUNK_SIZE})"
  )

  CompareNCFilesFamily (
    TEST_META_NAME ${TEST_BASE_NAME}_CHUNK_TYPE_vs_not_chunked
    FILE_META_NAME ${TEST_BASE_NAME}_output_CHUNK_TYPE.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.${RUN_T0}.nc
    MAGIC_STRING "CHUNK_TYPE"
    MAGIC_VALUES chunked not_chunked 
    LABELS rrtmgp physics
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_chunked_np${TEST_RANK_END}_omp1
                      ${FIXTURES_BASE_NAME}_not_chunked_np${TEST_RANK_END}_omp1)
endif()
