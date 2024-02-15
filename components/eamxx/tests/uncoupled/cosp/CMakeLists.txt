include (ScreamUtils)

# Create the test
CreateADUnitTest(cosp_standalone
  LABELS cosp physics
  LIBS eamxx_cosp
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
)

# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/USGS-gtopo30_ne4np4pg2_16x_converted.c20200527.nc)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)
