
# Set up a file to hold a list of all of the  tests
SET (HOMME_TEST_LIST ${CMAKE_BINARY_DIR}/tests/test_list.sh)
FILE (WRITE ${HOMME_TEST_LIST} "${POUND}!/bin/bash\n")
FILE (APPEND ${HOMME_TEST_LIST} "${POUND} A list of all of the HOMME tests\n")
FILE (APPEND ${HOMME_TEST_LIST} "\n")
SET (NUM_TEST_FILES 0)

SET(Homme_Build_DIR ${CMAKE_BINARY_DIR})
SET(Homme_Results_DIR ${CMAKE_SOURCE_DIR}/test/reg_test/results)

# Want this to be resetable...
IF (DEFINED HOMME_BASELINE_DIR)
  # Need to check that the files are there...
  SET(Homme_Baseline_DIR ${HOMME_BASELINE_DIR}/tests/baseline)
ELSE ()
  SET(Homme_Baseline_DIR ${CMAKE_BINARY_DIR}/tests/baseline)
ENDIF ()

MESSAGE(STATUS "Setting HOMME_BASELINE_DIR as ${Homme_Baseline_DIR}")

# Host specific
IF (Homme_Hostname STREQUAL "Yellowstone") 

  # Modify the results directory to point to the yellowstone results
  SET(Homme_Results_DIR ${Homme_Results_DIR}/yellowstone)

ELSEIF (Homme_OS STREQUAL "Darwin") 

  SET(Homme_Results_DIR ${Homme_Results_DIR}/darwin)

ENDIF ()

IF (NOT ${Homme_Submission_Type} STREQUAL none)
  IF (DEFINED HOMME_PROJID)
    MESSAGE(STATUS "Using HOMME Project ID (HOMME_PROJID) ${HOMME_PROJID}")
  ELSE ()
    # Check environment variables
    IF (DEFINED ENV{HOMME_PROJID})
      SET(HOMME_PROJID $ENV{HOMME_PROJID})
      MESSAGE(STATUS "Using HOMME Project ID (HOMME_PROJID) ${HOMME_PROJID} from environment")
    ELSE ()
      MESSAGE(FATAL_ERROR "No HOMME Project ID found. Please set HOMME_PROJID")
    ENDIF ()
  ENDIF ()
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
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/cmake_variables.sh
               ${CMAKE_BINARY_DIR}/tests/cmake_variables.sh @ONLY)

CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/testing-utils.sh
               ${CMAKE_BINARY_DIR}/tests/testing-utils.sh COPYONLY)
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/diffTol.py
               ${CMAKE_BINARY_DIR}/tests/diffTol.py COPYONLY)

# The first CMake test is running
#   this is a little awkward because submit_tests.sh uses 
#   submissions-list.sh which hasn't been created yet.
ADD_TEST(submitAndRunTests ${CMAKE_BINARY_DIR}/tests/submit_tests.sh "all")

ADD_CUSTOM_TARGET(check 
                  COMMAND ${CMAKE_CTEST_COMMAND} "--output-on-failure")

ADD_CUSTOM_TARGET(baseline
         COMMAND ${CMAKE_BINARY_DIR}/tests/submit_tests.sh "baseline")

IF (${HOMME_QUEUING})
  # The first CMake test is running
  #   this is a little awkward because submit_tests.sh uses 
  #   submissions-list.sh which hasn't been created yet.
  ADD_TEST(submitAndRunTests ${CMAKE_BINARY_DIR}/tests/submit_tests.sh "all")

  ADD_DEPENDENCIES(check submitAndRunTests)
  # Set the time limit to submit and run the tests at three hours
  SET_TESTS_PROPERTIES(submitAndRunTests PROPERTIES TIMEOUT 10800)

ENDIF ()

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

# End the file with the total number of test files
FILE (APPEND ${HOMME_TEST_LIST} "\n")
FILE (APPEND ${HOMME_TEST_LIST} "num_test_files=${NUM_TEST_FILES}\n")

# Now create the submission-list.sh
EXECUTE_PROCESS(
  COMMAND ${CMAKE_BINARY_DIR}/tests/createRunScripts.sh "all"
  RESULT_VARIABLE CREATE_RUNS_RESULT
  OUTPUT_VARIABLE CREATE_RUNS_OUTPUT
  ERROR_VARIABLE CREATE_RUNS_ERROR)

IF (NOT ${CREATE_RUNS_RESULT} EQUAL 0) 
  MESSAGE(FATAL_ERROR "The script createRunScripts.sh has failed with error
                        ${CREATE_RUNS_RESULT}")
ENDIF ()

