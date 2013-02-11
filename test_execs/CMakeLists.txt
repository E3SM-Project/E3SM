
# Set up a file to hold a list of all of the  tests
SET (HOMME_TEST_LIST ${CMAKE_BINARY_DIR}/tests/test_list.sh)
FILE (WRITE ${HOMME_TEST_LIST} "${POUND}!/bin/bash\n")
FILE (APPEND ${HOMME_TEST_LIST} "${POUND} A list of all of the HOMME tests\n")
FILE (APPEND ${HOMME_TEST_LIST} "\n")
SET (NUM_TEST_FILES 0)

IF(BUILD_HOMME_PREQX)
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

IF(BUILD_HOMME_SWEQX)
  ADD_SUBDIRECTORY(swtcA)
  ADD_SUBDIRECTORY(swtcB)
  ADD_SUBDIRECTORY(swtcC)
  createTest(swtc1)
  createTest(swtc5)
  createTest(swtc6)
ENDIF()

IF(BUILD_HOMME_SWDGX)
  ADD_SUBDIRECTORY(swtc-dgA)
  createTest(swtc1-dg)
  createTest(swtc2-dg)
  createTest(swtc5-dg)
ENDIF()

FILE (APPEND ${HOMME_TEST_LIST} "\n")
FILE (APPEND ${HOMME_TEST_LIST} "num_test_files=${NUM_TEST_FILES}\n")

SET(Homme_Build_DIR ${CMAKE_BINARY_DIR})
SET(Homme_Results_DIR ${CMAKE_SOURCE_DIR}/test/reg_test/results)

# Host specific
IF (Homme_Hostname STREQUAL "Yellowstone") 

  # Modify the results directory to point to the yellowstone results
  # TODO: deal with different compilers
  SET(Homme_Results_DIR ${CMAKE_SOURCE_DIR}/test/reg_test/results/yellowstone)
  
  CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/createRunScripts.sh 
                 ${CMAKE_BINARY_DIR}/tests/createRunScripts.sh @ONLY)
  CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/submit_lsf.sh
                 ${CMAKE_BINARY_DIR}/tests/submit_lsf.sh @ONLY)
  CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/run_tests/testing-utils.sh
                 ${CMAKE_BINARY_DIR}/tests/testing-utils.sh COPYONLY)
ENDIF ()
