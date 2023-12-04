INCLUDE (ScreamUtils)

# Create the test
CreateUnitTest(surface_coupling surface_coupling.cpp
  LABELS surface_coupling driver
  LIBS scream_control scream_share diagnostics
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

## Copy (and configure) yaml files needed by tests
set (RUN_T0 2021-10-12-45000)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)
