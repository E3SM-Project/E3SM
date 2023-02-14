#include(ScreamUtils)
include(ScreamUtils)

# The restart test requires cprnc
include (BuildCprnc)
BuildCprnc()

## Test output

configure_file(io_test_multisnap.yaml io_test_multisnap.yaml)
configure_file(io_test_template.yaml io_test_template.yaml)
configure_file(io_test_control.yaml io_test_control.yaml)

CreateUnitTest(io_test "io.cpp" "scream_io;diagnostics" LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

# Test output on SE grid
configure_file(io_test_se_grid.yaml io_test_se_grid.yaml)
CreateUnitTest(io_test_se_grid "io_se_grid.cpp" scream_io LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test output restart
# NOTE: Each restart test is a "setup" for the restart_check test,
# and cannot run in parallel with other restart tests,
# due to contention of the rpointer file

configure_file(io_test_restart.yaml io_test_restart.yaml)
configure_file(io_test_restart_check.yaml io_test_restart_check.yaml)
CreateUnitTest(output_restart_test "output_restart.cpp" scream_io LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  PROPERTIES RESOURCE_LOCK rpointer_file FIXTURES_SETUP restart_setup
)

foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
  set (SRC_FILE io_output_restart.AVERAGE.nsteps_x10.np${MPI_RANKS}.2000-01-01-00020.nc)
  set (TGT_FILE io_output_restart_check.AVERAGE.nsteps_x10.np${MPI_RANKS}.2000-01-01-00020.nc)
  add_test (NAME io_test_restart_check_np${MPI_RANKS}
            COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_property(TEST io_test_restart_check_np${MPI_RANKS}
               PROPERTY FIXTURES_REQUIRED restart_setup)
endforeach()

## Test remap output
CreateUnitTest(io_remap_test "io_remap_test.cpp" "scream_io;diagnostics" LABELS "io,remap"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

