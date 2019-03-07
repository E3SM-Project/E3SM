INCLUDE (ScreamUtils)

SET (num_mpi_ranks 1)
SET (config_defines "")

# Test dynamics folder
CreateUnitTest(dynamics "dynamics_tests.cpp" "scream_share" num_mpi_ranks config_defines)
