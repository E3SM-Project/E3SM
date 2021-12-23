#include(ScreamUtils)
include(ScreamUtils)

## Test output
CreateUnitTest(io_test "io.cpp" scream_io LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test restart
# Each restart test is a "setup" for the restart_check test,
# and cannot run in parallel with other restart tests
CreateUnitTest(restart_test "restart.cpp" scream_io LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  PROPERTIES RESOURCE_LOCK rpointer_file FIXTURES_SETUP restart_setup
)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/CmpNcdumpOutput.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/CmpNcdumpOutput.cmake @ONLY)

add_test (NAME io_test_restart_check
          COMMAND ${CMAKE_COMMAND} -P CmpNcdumpOutput.cmake
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  configure_file(io_test_instant.yaml io_test_instant_np${MPI_RANKS}.yaml)
  configure_file(io_test_average.yaml io_test_average_np${MPI_RANKS}.yaml)
  configure_file(io_test_max.yaml io_test_max_np${MPI_RANKS}.yaml)
  configure_file(io_test_min.yaml io_test_min_np${MPI_RANKS}.yaml)
  configure_file(io_test_multisnap.yaml io_test_multisnap_np${MPI_RANKS}.yaml)
  configure_file(io_test_restart.yaml io_test_restart_np${MPI_RANKS}.yaml)
endforeach()

# Make the restart check test run *after* the individual restart tests,
# as well as make it trigger them automatically
set_property(TEST io_test_restart_check PROPERTY FIXTURES_REQUIRED restart_setup)
