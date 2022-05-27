INCLUDE (ScreamUtils)

# Test atmosphere processes
if (NOT SCREAM_BASELINES_ONLY)

  SET (TEST_LABELS "rrtmgp;physics;driver")
  # Required libraries
  find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATH}/lib)
  set (NEED_LIBS scream_rrtmgp rrtmgp ${NETCDF_C} scream_control scream_share physics_share yakl)

  CreateUnitTest(
      rrtmgp_standalone_unit rrtmgp_standalone_unit.cpp "${NEED_LIBS}" LABELS ${TEST_LABELS}
      EXE_ARGS "--ekat-test-params rrtmgp_inputfile=${SCREAM_DATA_DIR}/init/rrtmgp-allsky.nc,rrtmgp_baseline=${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc"
  )
  # This test needs the allsky baselines file
  add_dependencies (rrtmgp_standalone_unit rrtmgp_allsky_baseline.nc)

  ## Create free running rrtmgp stand alone test that runs from an initial condition file.
  CreateUnitTest(
      rrtmgp_standalone "rrtmgp_standalone.cpp" "${NEED_LIBS}" LABELS ${TEST_LABELS}
      MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
      PROPERTIES FIXTURES_SETUP rrtmgp_generate_output_nc_files
  )

  # Copy yaml input file to run directory
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_unit.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_unit.yaml COPYONLY)

  # Ensure test input files are present in the data dir
  GetInputFile(init/screami_ne2np4L72_20220524.nc)

  # The RRTMGP stand-alone test that runs multi-step
  # Set AD configurable options
  SetVarDependingOnTestProfile(NUM_STEPS 2 5 48)
  set (ATM_TIME_STEP 1800)

  ## Copy (and configure) yaml files needed by tests
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
  configure_file(rrtmgp_standalone_output.yaml rrtmgp_standalone_output.yaml)

  ## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
  ## Only if max mpi ranks is >1
  # This test requires CPRNC
  if (SCREAM_TEST_MAX_RANKS GREATER 1)
    include (BuildCprnc)
    BuildCprnc()
    SET (BASE_TEST_NAME "rrtmgp")
    foreach (MPI_RANKS RANGE 2 ${SCREAM_TEST_MAX_RANKS})
      set (SRC_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.Steps_x${NUM_STEPS}.np${MPI_RANKS}.nc")
      set (TGT_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.Steps_x${NUM_STEPS}.np1.nc")
      set (TEST_NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb")
      add_test (NAME ${TEST_NAME}
                COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
                RESOURCE_LOCK ${BASE_TEST_NAME}
                FIXTURES_REQUIRED rrtmgp_generate_output_nc_files)
    endforeach()
  endif()
endif()
