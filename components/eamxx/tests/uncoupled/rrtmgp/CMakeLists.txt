INCLUDE (ScreamUtils)

# Test atmosphere processes
if (NOT SCREAM_BASELINES_ONLY)


  CreateUnitTest(rrtmgp_standalone_unit rrtmgp_standalone_unit.cpp
    LABELS rrtmgp physics driver
    LIBS scream_rrtmgp rrtmgp scream_control yakl diagnostics
    EXE_ARGS "--ekat-test-params rrtmgp_inputfile=${SCREAM_DATA_DIR}/init/rrtmgp-allsky.nc,rrtmgp_baseline=${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc"
  )
  # This test needs the allsky baselines file
  add_dependencies (rrtmgp_standalone_unit rrtmgp_allsky_baseline.nc)

  ## Create free running rrtmgp stand alone test that runs from an initial condition file.
  CreateUnitTest(rrtmgp_standalone "rrtmgp_standalone.cpp"
    LABELS rrtmgp physics driver
    LIBS scream_rrtmgp rrtmgp scream_control yakl diagnostics
    MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
    EXE_ARGS "--ekat-test-params inputfile=input.yaml"
    FIXTURES_SETUP_INDIVIDUAL rrtmgp_generate_output_nc_files
  )

  # Copy yaml input file to run directory
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_unit.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_unit.yaml)

  # Ensure test input files are present in the data dir
  GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})

  # The RRTMGP stand-alone test that runs multi-step
  # Set AD configurable options
  SetVarDependingOnTestSize(NUM_STEPS 2 5 48)
  set (ATM_TIME_STEP 1800)
  set (RUN_T0 2021-10-12-45000)

  ## Copy (and configure) yaml files needed by tests
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                  ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)
  set (COL_CHUNK_SIZE 1000)
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                  ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)

  ## Add a standalone test with chunked columns, and compare against non-chunked
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
  # Compute the chunk size, given number or ranks
  CreateUnitTestFromExec(
      rrtmgp_standalone_chunked rrtmgp_standalone
      LABELS rrtmgp physics driver
      MPI_RANKS ${TEST_RANK_END}
      EXE_ARGS "--ekat-test-params inputfile=input_chunked.yaml"
      FIXTURES_SETUP_INDIVIDUAL rrtmgp_chunked_generate_output
      PROPERTIES PASS_REGULAR_EXPRESSION "(beg.end: 0, ${COL_CHUNK_SIZE})"
  )

  # Compare chunked vs non-chunked radiation
  include (BuildCprnc)
  BuildCprnc()
  set (SRC_FILE "rrtmgp_standalone_output_chunked.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.${RUN_T0}.nc")
  set (TGT_FILE "rrtmgp_standalone_output.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.${RUN_T0}.nc")
  set (TEST_NAME "rrtmgp_chunked_vs_monolithic_bfb")
  add_test (NAME ${TEST_NAME}
            COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TEST_NAME} PROPERTIES
    LABELS "rrtmgp;physics;driver"
    FIXTURES_REQUIRED "rrtmgp_generate_output_nc_files_np${TEST_RANK_END}_omp1;rrtmgp_chunked_generate_output_np${TEST_RANK_END}_omp1")

  ## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
  ## Only if running with 2+ ranks configurations
  # This test requires CPRNC
  if (TEST_RANK_END GREATER TEST_RANK_START)
    math (EXPR CMP_RANK_START ${TEST_RANK_START}+1)
    set (FIXTURES_BASE_NAME rrtmgp_generate_output_nc_files)
    foreach (MPI_RANKS RANGE ${CMP_RANK_START} ${TEST_RANK_END})

      set (SRC_FILE "rrtmgp_standalone_output.INSTANT.nsteps_x${NUM_STEPS}.np${MPI_RANKS}.${RUN_T0}.nc")
      set (TGT_FILE "rrtmgp_standalone_output.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_START}.${RUN_T0}.nc")
      set (TEST_NAME "rrtmgp_np${TEST_RANK_START}_vs_np${MPI_RANKS}_bfb")
      add_test (NAME ${TEST_NAME}
                COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      set_tests_properties(${TEST_NAME} PROPERTIES
          LABELS "rrtmgp;physics;driver"
          FIXTURES_REQUIRED "${FIXTURES_BASE_NAME}_np${MPI_RANKS}_omp1;${FIXTURES_BASE_NAME}_np1_omp1")
    endforeach()
  endif()
endif()
