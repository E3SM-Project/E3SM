
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

# Determine if there is a submission system
OPTION(USE_QUEUING "Whether run tests through the queue" TRUE)
IF (USE_QUEUING)
  # lsf uses bsub
  FIND_PROGRAM(BSUB_PATH bsub) # Note this searches the environment
  IF (BSUB_PATH)
    MESSAGE(STATUS "Using LSF queuing")
    SET(HOMME_SUBMISSION_TYPE lsf)
    SET(HOMME_SUBMISSION_HEADER ${CMAKE_SOURCE_DIR}/cmake/queueSubHeaders/lsf.sh)
  ENDIF ()
  IF (NOT DEFINED HOMME_SUBMISSION_TYPE)
    # pbs uses qsub
    FIND_PROGRAM(QSUB_PATH qsub) # Note this searches the environment
    IF (QSUB_PATH)
      MESSAGE(STATUS "Using PBS queuing")
      SET(HOMME_SUBMISSION_TYPE pbs)
      SET(HOMME_SUBMISSION_HEADER ${CMAKE_SOURCE_DIR}/cmake/queueSubHeaders/pbs.sh)
    ENDIF ()
  ENDIF ()
  IF (NOT DEFINED HOMME_SUBMISSION_TYPE)
    # moab uses msub
    FIND_PROGRAM(MSUB_PATH msub) # Note this searches the environment
    IF (MSUB_PATH)
      MESSAGE(STATUS "Using MOAB queuing")
      SET(HOMME_SUBMISSION_TYPE moab)
      SET(HOMME_SUBMISSION_HEADER ${CMAKE_SOURCE_DIR}/cmake/queueSubHeaders/moab.sh)
    ENDIF ()
  ENDIF ()
  IF (NOT DEFINED HOMME_SUBMISSION_TYPE)
    # slurm uses squeue
    FIND_PROGRAM(SQUEUE_PATH squeue) # Note this searches the environment
    IF (SQUEUE_PATH)
      MESSAGE(STATUS "Using SLURM queuing")
      SET(HOMME_SUBMISSION_TYPE slurm)
      SET(HOMME_SUBMISSION_HEADER ${CMAKE_SOURCE_DIR}/cmake/queueSubHeaders/slurm.sh)
    ENDIF ()
  ENDIF ()
  IF (NOT DEFINED HOMME_SUBMISSION_TYPE)
    MESSAGE(STATUS "No queuing system found")
    SET(HOMME_SUBMISSION_TYPE none)
  ENDIF ()
ELSE ()
  MESSAGE(STATUS "Not searching for queuing system")
  SET(HOMME_SUBMISSION_TYPE none)
ENDIF ()

# Queuing specific info
IF (${HOMME_SUBMISSION_TYPE} STREQUAL none) 
  SET(HOMME_QUEUING FALSE)
  IF (USE_NUM_PROCS)
    MESSAGE(STATUS "Using ${USE_NUM_PROCS} processors")
  ELSE ()
    INCLUDE(ProcessorCount)
    ProcessorCount(MAX_NUM_PROCS)
    MESSAGE(STATUS "Warning: No queuing system. Maximum number of processes is ${MAX_NUM_PROCS}")
  ENDIF ()
ELSE ()
  IF (USE_NUM_PROCS)
    MESSAGE(STATUS "Warning: Variable USE_NUM_PROCS not supported when using queuing")
  ENDIF ()
  IF (DEFINED HOMME_PROJID)
    MESSAGE(STATUS "Using HOMME Project ID (HOMME_PROJID) ${HOMME_PROJID}")
  ELSE ()
    # Check environment variables
    IF (DEFINED ENV{HOMME_PROJID})
      SET(HOMME_PROJID $ENV{HOMME_PROJID})
      MESSAGE(STATUS "Using HOMME Project ID (HOMME_PROJID) ${HOMME_PROJID} from environment")
    ELSE ()
      MESSAGE(FATAL_ERROR "No HOMME Project ID found. Please set HOMME_PROJID
                           using -DHOMME_PROJID=XXXX0000 or in the environment")
    ENDIF ()
  ENDIF ()
  # We have queuing
  SET(HOMME_QUEUING TRUE)

  # Allow use to override the HOMME_SUBMISSION_HEADER
  IF (DEFINED SUBMISSION_HEADER)
    IF (EXISTS ${SUBMISSION_HEADER})
      SET(HOMME_SUBMISSION_HEADER ${SUBMISSION_HEADER})
    ELSE ()
      MESSAGE(FATAL_ERROR "SUBMISSION_HEADER set but file not found at: ${SUBMISSION_HEADER}")
    ENDIF ()
  ENDIF ()
  IF (DEFINED HOMME_SUBMISSION_HEADER)
    MESSAGE(STATUS "Using submission header ${HOMME_SUBMISSION_HEADER}")
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
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/check_baseline.sh
               ${CMAKE_BINARY_DIR}/tests/check_baseline.sh @ONLY)

CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/testing-utils.sh
               ${CMAKE_BINARY_DIR}/tests/testing-utils.sh COPYONLY)
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/diffTol.py
               ${CMAKE_BINARY_DIR}/tests/diffTol.py COPYONLY)

ADD_CUSTOM_TARGET(check 
                  COMMAND ${CMAKE_CTEST_COMMAND} "--output-on-failure")

# Force cprnc to be built when make check is run
ADD_DEPENDENCIES(check cprnc)

# Create a target for making the reference data
ADD_CUSTOM_TARGET(baseline
         COMMAND ${CMAKE_BINARY_DIR}/tests/submit_tests.sh "baseline")

# Create a target for making the reference data
ADD_TEST(verifyBaselineResults ${CMAKE_BINARY_DIR}/tests/check_baseline.sh)
ADD_DEPENDENCIES(check verifyBaselineResults)

IF (${HOMME_QUEUING})
  # The first CMake test is running
  #   this is a little awkward because submit_tests.sh uses 
  #   submissions-list.sh which hasn't been created yet.
  ADD_TEST(submitAndRunTests ${CMAKE_BINARY_DIR}/tests/submit_tests.sh "all")

  ADD_DEPENDENCIES(check submitAndRunTests)
  # Set the time limit to submit and run the tests at three hours
  SET_TESTS_PROPERTIES(submitAndRunTests PROPERTIES TIMEOUT 10800)

ENDIF ()

# Add the test exec subdirs for the sweqx executable
IF(${BUILD_HOMME_SWEQX})
  ADD_SUBDIRECTORY(swtcA)
  ADD_SUBDIRECTORY(swtcB)
ENDIF()


# Add the test exec subdirs for the preqx executable
IF(${BUILD_HOMME_PREQX})
  ADD_SUBDIRECTORY(baroC)
  ADD_SUBDIRECTORY(baroCam)
  ADD_SUBDIRECTORY(preqx-nlev20-interp)
  ADD_SUBDIRECTORY(preqx-nlev20-native)
  ADD_SUBDIRECTORY(preqx-nlev30-interp)
  ADD_SUBDIRECTORY(preqx-nlev30-native)
  ADD_SUBDIRECTORY(preqx-nlev60-interp)
  ADD_SUBDIRECTORY(preqx-nlev60-native)
  IF(${BUILD_HOMME_PREQX_ACC})
    ADD_SUBDIRECTORY(baroCam-acc)
  ENDIF()
ENDIF()

# Add the test exec subdirs for the prim executable
IF(${BUILD_HOMME_PESE})
  ADD_SUBDIRECTORY(pese-nlev20)
  ADD_SUBDIRECTORY(pese-nlev30)
  ADD_SUBDIRECTORY(pese-nlev60)
ENDIF()

IF(${BUILD_HOMME_THETA})
  ADD_SUBDIRECTORY(theta-nlev20)
  ADD_SUBDIRECTORY(theta-nlev30)
  ADD_SUBDIRECTORY(theta-nlev60)
ENDIF()

# Add the test exec subdirs for the prim executable
IF(${BUILD_HOMME_PRIM})
  ADD_SUBDIRECTORY(primA)
ENDIF()

# Add the test exec subdirs for the swim executable
IF(${BUILD_HOMME_SWIM})
  ADD_SUBDIRECTORY(swim5)
ENDIF()

# Read the test-list.cmake file to get the HOMME_TESTS list
INCLUDE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/test-list.cmake)
# Need a properly scoped variable to store the test directories for the
#  replaceBaselineResults target
SET(TEST_DIR_LIST "")
createTests(HOMME_TESTS)

IF(${BUILD_HOMME_PREQX} AND ${BUILD_HOMME_PREQX_ACC})
  # add a test dependency here:
  set_property(TEST baroCamMoist-acc APPEND PROPERTY DEPENDS baroCamMoist)
ENDIF()

# End the file with the total number of test files
FILE (APPEND ${HOMME_TEST_LIST} "\n")
FILE (APPEND ${HOMME_TEST_LIST} "NUM_TEST_FILES=${NUM_TEST_FILES}\n")

# Create a target for making the reference data
ADD_CUSTOM_TARGET(replaceBaselineResults)
FOREACH(TEST_DIR ${TEST_DIR_LIST})
  GET_FILENAME_COMPONENT(TEST_BASE_NAME ${TEST_DIR} NAME)
  ADD_CUSTOM_COMMAND(TARGET replaceBaselineResults PRE_BUILD
           COMMAND ${CMAKE_COMMAND} -E copy_directory 
           ${TEST_DIR} ${Homme_Baseline_DIR}/${TEST_BASE_NAME})
  ADD_CUSTOM_COMMAND(TARGET replaceBaselineResults PRE_BUILD
           COMMAND ${CMAKE_COMMAND} -E echo 
           "Copying test results ${TEST_DIR} to baseline directory ${Homme_Baseline_DIR}/${TEST_BASE_NAME}")
ENDFOREACH()

# Now create the submission-list.sh
EXECUTE_PROCESS(
  COMMAND ${CMAKE_BINARY_DIR}/tests/createRunScripts.sh "all"
  RESULT_VARIABLE CREATE_RUNS_RESULT
  OUTPUT_VARIABLE CREATE_RUNS_OUTPUT
  ERROR_VARIABLE CREATE_RUNS_ERROR)

IF (NOT ${CREATE_RUNS_RESULT} EQUAL 0) 
  MESSAGE(FATAL_ERROR "The script \"createRunScripts.sh all\" has failed with error
                        ${CREATE_RUNS_RESULT}\n"
          "OUTPUT: \n ${CREATE_RUNS_OUTPUT}\n"
          "ERROR: \n ${CREATE_RUNS_ERROR}")
ENDIF ()

EXECUTE_PROCESS(
  COMMAND ${CMAKE_BINARY_DIR}/tests/createRunScripts.sh "baseline"
  RESULT_VARIABLE CREATE_RUNS_RESULT
  OUTPUT_VARIABLE CREATE_RUNS_OUTPUT
  ERROR_VARIABLE CREATE_RUNS_ERROR)

IF (NOT ${CREATE_RUNS_RESULT} EQUAL 0) 
  MESSAGE(FATAL_ERROR "The script \"createRunScripts.sh baseline\" has failed with error
                        ${CREATE_RUNS_RESULT}\n"
          "OUTPUT: \n ${CREATE_RUNS_OUTPUT}\n"
          "ERROR: \n ${CREATE_RUNS_ERROR}")

ENDIF ()

