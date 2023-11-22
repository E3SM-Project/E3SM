include (ScreamUtils)

# Create the test
CreateADUnitTest(mam4_nucleation_standalone
  LABELS mam4 physics
  LIBS mam
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
)

# Set AD configurable options
set (ATM_TIME_STEP 1)
SetVarDependingOnTestSize(NUM_STEPS 12 24 36)
set (RUN_T0 2021-10-12-45000)

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
