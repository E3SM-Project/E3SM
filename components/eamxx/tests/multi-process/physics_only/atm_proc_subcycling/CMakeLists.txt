INCLUDE (ScreamUtils)

# Create the exec
CreateADUnitTestExec (shoc_p3
  LIBS shoc p3 diagnostics)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

set (RUN_T0 2021-10-12-45000)

# Run a test with subcycling
set (NUM_SUBCYCLES 3)
set (NUM_STEPS  1)
set (POSTFIX subcycled)
set (ATM_TIME_STEP 900)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_subcycled.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_subcycled.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_tend.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_tend_subcycled.yaml)

CreateUnitTestFromExec (shoc_p3_subcycled shoc_p3
      EXE_ARGS "--args -ifile=input_subcycled.yaml"
      FIXTURES_SETUP shoc_p3_subcycled)

# Run a test without subcycling and more steps
set (NUM_SUBCYCLES 1)
set (NUM_STEPS  3)
set (ATM_TIME_STEP 300)
set (POSTFIX monolithic)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_monolithic.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_monolithic.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_tend.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_tend_monolithic.yaml)
CreateUnitTestFromExec (shoc_p3_monolithic shoc_p3
      EXE_ARGS "--args -ifile=input_monolithic.yaml"
      FIXTURES_SETUP shoc_p3_monolithic)

# Finally, compare output of the two tests
include (BuildCprnc)
BuildCprnc()

set (SRC_FILE "shoc_p3_monolithic.INSTANT.nsteps_x3.np1.${RUN_T0}.nc")
set (TGT_FILE "shoc_p3_subcycled.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
set (TEST_NAME check_subcycling)
add_test (NAME ${TEST_NAME}
          COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(${TEST_NAME} PROPERTIES
      LABELS shoc p3 infrastructure
      FIXTURES_REQUIRED "shoc_p3_monolithic;shoc_p3_subcycled")

# Check calculation of shoc tendencies when the parent group is subcycled
set (script ${SCREAM_BASE_DIR}/scripts/check-tendencies)
set (fname shoc_p3_tend_subcycled.INSTANT.nsteps_x1.np1.${RUN_T0}.nc)
add_test (NAME ${TEST_NAME}_tend_check
          COMMAND ${script} -f ${fname} -v tke -t shoc_tke_tend
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (${TEST_NAME}_tend_check PROPERTIES
      LABELS shoc p3 infrastructure
      FIXTURES_REQUIRED shoc_p3_subcycled)
