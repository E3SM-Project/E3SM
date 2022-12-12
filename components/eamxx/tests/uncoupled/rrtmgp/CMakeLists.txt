INCLUDE (ScreamUtils)

# Test atmosphere processes
if (NOT SCREAM_BASELINES_ONLY)

  SET (TEST_LABELS "rrtmgp;physics;driver")
  # Required libraries
  find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATH}/lib)
  set (NEED_LIBS scream_rrtmgp rrtmgp ${NETCDF_C} scream_control scream_share physics_share yakl diagnostics)

  CreateUnitTest(
      rrtmgp_standalone_unit rrtmgp_standalone_unit.cpp "${NEED_LIBS};rrtmgp_test_utils" LABELS ${TEST_LABELS}
      EXE_ARGS "--ekat-test-params rrtmgp_inputfile=${SCREAM_DATA_DIR}/init/rrtmgp-allsky.nc,rrtmgp_baseline=${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc"
  )
  # This test needs the allsky baselines file
  add_dependencies (rrtmgp_standalone_unit rrtmgp_allsky_baseline.nc)

  ## Create free running rrtmgp stand alone test that runs from an initial condition file.
  CreateUnitTest(
      rrtmgp_standalone "rrtmgp_standalone.cpp" "${NEED_LIBS}" LABELS ${TEST_LABELS}
      MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
      EXE_ARGS "--ekat-test-params inputfile=input.yaml"
      PROPERTIES FIXTURES_SETUP rrtmgp_generate_output_nc_files
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

  ## Copy (and configure) yaml files needed by tests
  configure_file (rrtmgp_standalone_output.yaml rrtmgp_standalone_output.yaml)
  set (COL_CHUNK_SIZE 1000)
  configure_file (input.yaml input.yaml)

  ## Add a standalone test with chunked columns, and compare against non-chunked
  set (SUFFIX "_chunked")
  math (EXPR COL_PER_RANK "218 / ${TEST_RANK_END}")
  math (EXPR COL_CHUNK_SIZE "${COL_PER_RANK} / 2")
  if (COL_CHUNK_SIZE LESS 1)
    message (FATAL_ERROR "Error! Chunk size for rrtmgp unit test is less than 1.")
  endif()

  configure_file (input.yaml input_chunked.yaml)
  configure_file (rrtmgp_standalone_output.yaml rrtmgp_standalone_output_chunked.yaml)
  # Compute the chunk size, given number or ranks
  CreateUnitTestFromExec(
      rrtmgp_standalone_chunked rrtmgp_standalone
      LABELS ${TEST_LABELS}
      MPI_RANKS ${TEST_RANK_END}
      EXE_ARGS "--ekat-test-params inputfile=input_chunked.yaml"
      PROPERTIES FIXTURES_SETUP rrtmgp_chunked_generate_output
      PASS_REGULAR_EXPRESSION "(beg.end: 0, ${COL_CHUNK_SIZE})"
  )

  # Compare chunked vs non-chunked radiation
  include (BuildCprnc)
  BuildCprnc()
  set (SRC_FILE "rrtmgp_standalone_output_chunked.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.nc")
  set (TGT_FILE "rrtmgp_standalone_output.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.nc")
  set (TEST_NAME "rrtmgp_chunked_vs_monolithic_bfb")
  add_test (NAME ${TEST_NAME}
            COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
            FIXTURES_REQUIRED "rrtmgp_generate_output_nc_files;rrtmgp_chunked_generate_output")

  ## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
  ## Only if running with 2+ ranks configurations
  # This test requires CPRNC
  if (TEST_RANK_END GREATER TEST_RANK_START)
    foreach (MPI_RANKS RANGE ${TEST_RANK_START} ${TEST_RANK_END})
      set (SRC_FILE "rrtmgp_standalone_output.INSTANT.nsteps_x${NUM_STEPS}.np${MPI_RANKS}.nc")
      set (TGT_FILE "rrtmgp_standalone_output.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_START}.nc")
      set (TEST_NAME "rrtmgp_np${TEST_RANK_START}_vs_np${MPI_RANKS}_bfb")
      add_test (NAME ${TEST_NAME}
                COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
                FIXTURES_REQUIRED rrtmgp_generate_output_nc_files)
    endforeach()
  endif()
endif()
