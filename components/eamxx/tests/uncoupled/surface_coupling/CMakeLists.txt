INCLUDE (ScreamUtils)

# Create the test
set (TEST_LABELS "surface_coupling;driver")
set (NEED_LIBS scream_control scream_share diagnostics)
CreateUnitTest(surface_coupling surface_coupling.cpp "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  PROPERTIES FIXTURES_SETUP surface_coupling_generate_output_nc_files
)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/USGS-gtopo30_ne4np4_16x.c20160612.nc)

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/surface_coupling_output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/surface_coupling_output.yaml)
