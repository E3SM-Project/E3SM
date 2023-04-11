include (ScreamUtils)

# Create the test
SET (TEST_LABELS "spa;physics;driver")
set (NEED_LIBS spa scream_control scream_share diagnostics)

# Test atmosphere processes
CreateUnitTest(spa_stand_alone "spa_stand_alone.cpp" "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  PROPERTIES FIXTURES_SETUP spa_generate_output_nc_files
)

# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 2 4 48)
set (ATM_TIME_STEP 864000)

## Copy (and configure) yaml files needed by tests
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_standalone_output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/spa_standalone_output.yaml)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  ${EAMxx_tests_IC_FILE_72lev}
  map_ne4np4_to_ne2np4_mono.nc
  spa_file_unified_and_complete_ne4_20220428.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(scream/init/${file})
endforeach()

## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if running with 2+ ranks configurations
# This test requires CPRNC
if (TEST_RANK_END GREATER TEST_RANK_START)
  include (BuildCprnc)
  BuildCprnc()
  SET (BASE_TEST_NAME "spa")
  foreach (MPI_RANKS RANGE ${TEST_RANK_START} ${TEST_RANK_END})
    set (SRC_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.nsteps_x1.np${MPI_RANKS}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.nsteps_x1.np${TEST_RANK_START}.nc")
    set (TEST_NAME "${BASE_TEST_NAME}_np${TEST_RANK_START}_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED spa_generate_output_nc_files)
  endforeach()
endif()
