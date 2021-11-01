include (ScreamUtils)

# Create the test
SET (TEST_LABELS "spa;physics")
set (NEED_LIBS spa scream_control scream_share)
CreateUnitTest(spa_stand_alone "spa_stand_alone.cpp" "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  PROPERTIES FIXTURES_SETUP spa_main
)

# Set AD configurable options
SetVarDependingOnTestProfile(NUM_STEPS 2 4 48)
set (ATM_TIME_STEP 864000)

## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_np${MPI_RANKS}.yaml)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_standalone_output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/spa_standalone_output_np${MPI_RANKS}.yaml)

endforeach()

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
    set (SRC_FILE "${BASE_TEST_NAME}_standalone_output_np${MPI_RANKS}.INSTANT.Steps_x${NUM_STEPS}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_standalone_output_np1.INSTANT.Steps_x${NUM_STEPS}.nc")
    configure_file (${SCREAM_BASE_DIR}/cmake/CprncTest.cmake
                    ${CMAKE_CURRENT_BINARY_DIR}/CprncTest_np${MPI_RANKS}.cmake @ONLY)
    set(TEST_NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P CprncTest_np${MPI_RANKS}.cmake
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED spa_main)
  endforeach()
endif()
