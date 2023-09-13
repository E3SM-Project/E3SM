INCLUDE (ScreamUtils)

# Create the exec
CreateUnitTestExec (shoc_p3_nudging shoc_p3_nudging.cpp "shoc;p3;nudging;scream_control;diagnostics")

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

set (RUN_T0 2021-10-12-45000)

# Run a test to setup nudging source data:
set (NUM_STEPS  5)
set (POSTFIX source_data)
set (ATM_TIME_STEP 300)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_source_data.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_source_data.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_source_data.yaml)
CreateUnitTestFromExec (shoc_p3_source shoc_p3_nudging
      EXE_ARGS "--use-colour no --ekat-test-params ifile=input_source_data.yaml"
      PROPERTIES FIXTURES_SETUP source_data)

# Run a test with nudging turned on:
set (NUM_STEPS  5)
set (ATM_TIME_STEP 300)
set (POSTFIX nudged)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_nudging.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_nudging.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_nudged.yaml)
CreateUnitTestFromExec (shoc_p3_nudged shoc_p3_nudging
      EXE_ARGS "--use-colour no --ekat-test-params ifile=input_nudging.yaml"
      PROPERTIES FIXTURES_SETUP nudged FIXTURES_REQUIRED source_data)

# Finally, compare output of the two tests
include (BuildCprnc)
BuildCprnc()

set (SRC_FILE "shoc_p3_source_data.INSTANT.nsteps_x3.${RUN_T0}.nc")
set (TGT_FILE "shoc_p3_nudged.INSTANT.nsteps_x1.${RUN_T0}.nc")
set (TEST_NAME check_subcycling)
add_test (NAME ${TEST_NAME}
          COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "shoc;p3"
          FIXTURES_REQUIRED "source_data;nudged")
