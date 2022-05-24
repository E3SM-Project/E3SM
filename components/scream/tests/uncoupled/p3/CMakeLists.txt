include (ScreamUtils)

# Create the test
SET (TEST_LABELS "p3;physics;driver")
set (NEED_LIBS p3 scream_control scream_share)
CreateUnitTest(p3_standalone "p3_standalone.cpp" "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  PROPERTIES FIXTURES_SETUP p3_generate_output_nc_files
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestProfile(NUM_STEPS 2 5 48)  # 1h 2.5h 24h

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(p3_standalone_output.yaml p3_standalone_output.yaml)

# Ensure test input files are present in the data dir
GetInputFile(init/screami_ne2np4L72_20220524.nc)

## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if max mpi ranks is >1
# This test requires CPRNC
if (SCREAM_TEST_MAX_RANKS GREATER 1)
  include (BuildCprnc)
  BuildCprnc()
  SET (BASE_TEST_NAME "p3")
  foreach (MPI_RANKS RANGE 2 ${SCREAM_TEST_MAX_RANKS})
    set (SRC_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.Steps_x${NUM_STEPS}.np${MPI_RANKS}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.Steps_x${NUM_STEPS}.np1.nc")
    set (TEST_NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED p3_generate_output_nc_files)
  endforeach()
endif()
