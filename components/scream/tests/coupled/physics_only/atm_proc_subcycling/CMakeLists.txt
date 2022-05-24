INCLUDE (ScreamUtils)

# Create the exec
CreateUnitTestExec (shoc_p3 shoc_p3.cpp "shoc;p3;scream_control")

# Ensure test input files are present in the data dir
GetInputFile(init/screami_ne2np4L72_20220524.nc)

# Run a test with subcycling
set (NUM_SUBCYCLES 3)
set (NUM_STEPS  1)
set (POSTFIX subcycled)
set (ATM_TIME_STEP 900)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_subcycled.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_subcycled.yaml)
CreateUnitTestFromExec (shoc_p3_subcycled shoc_p3
      EXE_ARGS "--use-colour no --ekat-test-params ifile=input_subcycled.yaml"
      PROPERTIES FIXTURES_SETUP subcycled)

# Run a test without subcycling and more steps
set (NUM_SUBCYCLES 1)
set (NUM_STEPS  3)
set (ATM_TIME_STEP 300)
set (POSTFIX monolithic)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_monolithic.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_monolithic.yaml)
CreateUnitTestFromExec (shoc_p3_monolithic shoc_p3
      EXE_ARGS "--use-colour no --ekat-test-params ifile=input_monolithic.yaml"
      PROPERTIES FIXTURES_SETUP monolithic)

# Finally, compare output of the two tests
include (BuildCprnc)
BuildCprnc()

set (SRC_FILE "shoc_p3_monolithic.INSTANT.Steps_x3.nc")
set (TGT_FILE "shoc_p3_subcycled.INSTANT.Steps_x1.nc")
set (TEST_NAME check_subcycling)
add_test (NAME ${TEST_NAME}
          COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(${TEST_NAME} PROPERTIES LABELS "shoc;p3"
          FIXTURES_REQUIRED "monolithic;subcycled")
