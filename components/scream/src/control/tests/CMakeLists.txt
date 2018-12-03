INCLUDE (ScreamUtils)

SET (num_mpi_ranks 1)
SET (config_defines "")

# Test control folder
CreateUnitTest(control "control_tests.cpp" "scream_share;scream_control" num_mpi_ranks config_defines)
