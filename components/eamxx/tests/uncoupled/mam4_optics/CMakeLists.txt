include (ScreamUtils)

# Create the test
SET (TEST_LABELS "mam4;physics;driver")
set (NEED_LIBS mam scream_control scream_share diagnostics)
CreateUnitTest(mam4_optics_standalone "mam4_optics_standalone.cpp" "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_START}
  PROPERTIES FIXTURES_SETUP mam4_generate_output_nc_files
)

# Set AD configurable options
set (ATM_TIME_STEP 1)
SetVarDependingOnTestSize(NUM_STEPS 12 24 36)
set (RUN_T0 2021-10-12-45000)

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(mam4_optics_output.yaml mam4_optics_output.yaml)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
