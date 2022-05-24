include (ScreamUtils)

# Create the test
SET (TEST_LABELS "shoc;physics;driver")
set (NEED_LIBS shoc scream_control scream_share)
CreateUnitTest(shoc_standalone "shoc_standalone.cpp" "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  PROPERTIES FIXTURES_SETUP shoc_generate_output_nc_files
)

# Set AD configurable options
SetVarDependingOnTestProfile(NUM_STEPS 2 5 48)
set (ATM_TIME_STEP 1800)

# Determine num subcycles needed to keep shoc dt<=300s
set (SHOC_MAX_DT 300)
math (EXPR NUM_SUBCYCLES "(${ATM_TIME_STEP} + ${SHOC_MAX_DT} - 1) / ${SHOC_MAX_DT}")

# Ensure test input files are present in the data dir
GetInputFile(init/screami_ne2np4L72_20220524.nc)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(shoc_standalone_output.yaml shoc_standalone_output.yaml)

## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if max mpi ranks is >1
# This test requires CPRNC
if (SCREAM_TEST_MAX_RANKS GREATER 1)
  include (BuildCprnc)
  BuildCprnc()
  SET (BASE_TEST_NAME "shoc")
  foreach (MPI_RANKS RANGE 2 ${SCREAM_TEST_MAX_RANKS})
    set (SRC_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.Steps_x${NUM_STEPS}.np${MPI_RANKS}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.Steps_x${NUM_STEPS}.np1.nc")
    set (TEST_NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME} 
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED shoc_generate_output_nc_files)
  endforeach()
endif()
