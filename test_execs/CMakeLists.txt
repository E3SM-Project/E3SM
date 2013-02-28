
# Set up a file to hold a list of all of the  tests
SET (HOMME_TEST_LIST ${CMAKE_BINARY_DIR}/tests/test_list.sh)
FILE (WRITE ${HOMME_TEST_LIST} "${POUND}!/bin/bash\n")
FILE (APPEND ${HOMME_TEST_LIST} "${POUND} A list of all of the HOMME tests\n")
FILE (APPEND ${HOMME_TEST_LIST} "\n")
SET (NUM_TEST_FILES 0)

SET(Homme_Build_DIR ${CMAKE_BINARY_DIR})
SET(Homme_Results_DIR ${CMAKE_SOURCE_DIR}/test/reg_test/results)

# Host specific
IF (Homme_Hostname STREQUAL "Yellowstone") 

  # Modify the results directory to point to the yellowstone results
  # TODO: deal with different compilers
  SET(Homme_Results_DIR ${Homme_Results_DIR}/yellowstone)

  SET(Homme_Submission_Type lsf)

ELSEIF (Homme_OS STREQUAL "Darwin") 

  SET(Homme_Results_DIR ${Homme_Results_DIR}/darwin)

  SET(Homme_Submission_Type std)

ENDIF ()

# Because of the queuing system on machines like Yellowstone the testing can be a little awkward
# The following scripts prepare submission scripts and allow submission to
# queues on machines that have queues.
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/createRunScripts.sh 
               ${CMAKE_BINARY_DIR}/tests/createRunScripts.sh @ONLY)
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/submit_tests.sh
               ${CMAKE_BINARY_DIR}/tests/submit_tests.sh @ONLY)
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/diff_output.sh
               ${CMAKE_BINARY_DIR}/tests/diff_output.sh @ONLY)
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/testing-utils.sh
               ${CMAKE_BINARY_DIR}/tests/testing-utils.sh COPYONLY)

# The first CMake test is running
#   this is a little awkward because submit_tests.sh uses 
#   submissions-list.sh which hasn't been created yet.
ADD_TEST(submitAndRunTests ${CMAKE_BINARY_DIR}/tests/submit_tests.sh)

# Set the time limit to submit and run the tests at three hours
SET_TESTS_PROPERTIES(submitAndRunTests PROPERTIES TIMEOUT 10800)

IF(${BUILD_HOMME_PREQX})
  ADD_SUBDIRECTORY(baroA)
  ADD_SUBDIRECTORY(baroB)
  ADD_SUBDIRECTORY(baroC)
  createTest(baro1a)
  createTest(baro1b)
  createTest(baro2a)
  createTest(baro2b)
  createTest(baro2c)
  createTest(baro2d)
ENDIF()

IF(${BUILD_HOMME_SWEQX})
  ADD_SUBDIRECTORY(swtcA)
  ADD_SUBDIRECTORY(swtcB)
  ADD_SUBDIRECTORY(swtcC)
  createTest(swtc1)
  createTest(swtc5)
  createTest(swtc6)
ENDIF()

IF(${BUILD_HOMME_SWDGX})
  ADD_SUBDIRECTORY(swtc-dgA)
  createTest(swtc1-dg)
  createTest(swtc2-dg)
  createTest(swtc5-dg)
ENDIF()

# SNow create the submission-list.sh
FILE (APPEND ${HOMME_TEST_LIST} "\n")
FILE (APPEND ${HOMME_TEST_LIST} "num_test_files=${NUM_TEST_FILES}\n")

# Now create the submission-list.sh
EXECUTE_PROCESS(
  COMMAND ${CMAKE_BINARY_DIR}/tests/createRunScripts.sh
  RESULT_VARIABLE CREATE_RUNS_RESULT
  OUTPUT_VARIABLE CREATE_RUNS_OUTPUT
  ERROR_VARIABLE CREATE_RUNS_ERROR
)


