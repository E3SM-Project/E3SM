#include(ScreamUtils)
include(ScreamUtils)

set (ARRAY_SCORPIO_SRCS
  scorpio_tests.cpp
  io.cpp
  restart.cpp
)

## Test atmosphere processes
# TODO:
# This commented test may be needed for future development.
# Currently tests a functionality that isn't yet in the code.
# Keeping for now for use when that development task is finished.
#CreateUnitTest(scorpio_interface "scorpio_tests.cpp" scream_io 
#  MPI_RANKS 1 4
#)
#configure_file(scorpio_output_baseline.nc scorpio_output_baseline.nc COPYONLY)

CreateUnitTest(io_test "io.cpp" scream_io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

CreateUnitTest(restart_test "restart.cpp" scream_io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

# As soon as changes from E3SM-Project/EKAT#79 are integrated in SCREAM,
# modify CreateUnitTest to make the following handled in the macro.
set(tests_names)
foreach(RANK RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  list(APPEND tests_names restart_test_ut_np${RANK}_omp1)
endforeach()
set_tests_properties (${tests_names} PROPERTIES RESOURCE_LOCK rpointer_file)

## Copy (and configure) yaml files needed by tests
foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  configure_file(io_test_instant.yaml io_test_instant_np${MPI_RANKS}.yaml)
  configure_file(io_test_average.yaml io_test_average_np${MPI_RANKS}.yaml)
  configure_file(io_test_max.yaml io_test_max_np${MPI_RANKS}.yaml)
  configure_file(io_test_min.yaml io_test_min_np${MPI_RANKS}.yaml)
  configure_file(io_test_restart.yaml io_test_restart_np${MPI_RANKS}.yaml)
endforeach()
