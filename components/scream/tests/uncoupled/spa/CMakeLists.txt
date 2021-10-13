include (ScreamUtils)

# Create the test
SET (TEST_LABELS "spa;physics")
set (NEED_LIBS spa scream_control scream_share)
# Make sure the SPA data directory is correctly defined
if (NOT EXISTS ${SCREAM_SPA_DATA_DIR})
  message(FATAL_ERROR "Failed, SCREAM_SPA_DATA_DIR = ${SCREAM_SPA_DATA_DIR} does not exist")
endif()
if (NOT IS_DIRECTORY ${SCREAM_SPA_DATA_DIR})
  message(FATAL_ERROR "Failed, SCREAM_SPA_DATA_DIR = ${SCREAM_SPA_DATA_DIR} is not a directory")
endif() 
if (NOT IS_ABSOLUTE ${SCREAM_SPA_DATA_DIR})
  message(FATAL_ERROR "Failed, SCREAM_SPA_DATA_DIR = ${SCREAM_SPA_DATA_DIR} is not an absolute directory")
endif() 
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

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/spa_init_ne2np4.nc COPYONLY)

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
    configure_file (${SCREAM_BASE_DIR}/cmake/CprncTest.cmake
                    ${CMAKE_CURRENT_BINARY_DIR}/CprncTest_np${MPI_RANKS}.cmake @ONLY)
    set(TEST_NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P CprncTest_np${MPI_RANKS}.cmake
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED spa_generate_output_nc_files)
  endforeach()
endif()
