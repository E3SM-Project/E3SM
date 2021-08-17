#include(ScreamUtils)
include(ScreamUtils)

set (ARRAY_SCORPIO_SRCS
  io.cpp
  restart.cpp
)

## Test atmosphere processes
CreateUnitTest(io_test "io.cpp" scream_io LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

CreateUnitTest(restart_test "restart.cpp" scream_io LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

# As soon as changes from E3SM-Project/EKAT#79 are integrated in SCREAM,
# modify CreateUnitTest to make the following handled in the macro.
set(tests_names)
foreach(RANK RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  list(APPEND tests_names restart_test_ut_np${RANK}_omp1)
endforeach()
set_tests_properties (${tests_names} PROPERTIES RESOURCE_LOCK rpointer_file)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/CmpNcdumpOutput.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/CmpNcdumpOutput.cmake @ONLY)

add_test (NAME io_test_restart_check
          COMMAND ${CMAKE_COMMAND} -P CmpNcdumpOutput.cmake
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  # Each restart test is a "setup" for the restart_check test, that runs cprnc on their outputs.
  set_property(TEST restart_test_ut_np${MPI_RANKS}_omp1 PROPERTY FIXTURES_SETUP restart_setup)

  configure_file(io_test_instant.yaml io_test_instant_np${MPI_RANKS}.yaml)
  configure_file(io_test_average.yaml io_test_average_np${MPI_RANKS}.yaml)
  configure_file(io_test_max.yaml io_test_max_np${MPI_RANKS}.yaml)
  configure_file(io_test_min.yaml io_test_min_np${MPI_RANKS}.yaml)
  configure_file(io_test_restart.yaml io_test_restart_np${MPI_RANKS}.yaml)
  configure_file(io_test_restart_check.yaml io_test_restart_check_np${MPI_RANKS}.yaml)

endforeach()

# Make the restart check test run *after* the individual restart tests,
# as well as make it trigger them automatically
set_property(TEST io_test_restart_check PROPERTY FIXTURES_REQUIRED restart_setup)
