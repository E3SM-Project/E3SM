INCLUDE (ScreamUtils)

SET (num_mpi_ranks 1)
SET (config_defines "")

# Test control folder
CreateUnitTest(ping_pong "ping_pong_test.cpp" "scream_control;scream_share" num_mpi_ranks config_defines)

# Unit test the ad
CreateUnitTest(ad_ut "ad_tests.cpp" "scream_control;scream_share" num_mpi_ranks config_defines)
