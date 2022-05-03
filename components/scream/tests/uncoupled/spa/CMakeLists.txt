include (ScreamUtils)

# Create the test
SET (TEST_LABELS "spa;physics")
set (NEED_LIBS spa scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(spa_stand_alone "spa_stand_alone.cpp" "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  PROPERTIES FIXTURES_SETUP spa_generate_output_nc_files
)

# Set AD configurable options
SetVarDependingOnTestProfile(NUM_STEPS 2 4 48)
set (ATM_TIME_STEP 864000)

## Copy (and configure) yaml files needed by tests
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_standalone_output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/spa_standalone_output.yaml)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  spa_init_ne2np4.nc
  map_ne4np4_to_ne2np4_mono.nc
  spa_file_unified_and_complete_ne4_20220428.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(init/${file})
endforeach()

## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if max mpi ranks is >1
# This test requires CPRNC
if (SCREAM_TEST_MAX_RANKS GREATER 1)
  include (BuildCprnc)
  BuildCprnc()
  SET (BASE_TEST_NAME "spa")
  foreach (MPI_RANKS RANGE 2 ${SCREAM_TEST_MAX_RANKS})
    set (SRC_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.Steps_x${NUM_STEPS}.np${MPI_RANKS}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.Steps_x${NUM_STEPS}.np1.nc")
    set (TEST_NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED spa_generate_output_nc_files)
  endforeach()
endif()
