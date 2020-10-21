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

CreateUnitTest(io_test "io.cpp" scream_io)

CreateUnitTest(restart_test "restart.cpp" scream_io)

## Copy files needed by tests
configure_file(io_test_instant.yaml io_test_instant.yaml COPYONLY)
configure_file(io_test_average.yaml io_test_average.yaml COPYONLY)
configure_file(io_test_max.yaml io_test_max.yaml COPYONLY)
configure_file(io_test_min.yaml io_test_min.yaml COPYONLY)
configure_file(io_test_restart.yaml io_test_restart.yaml COPYONLY)
