INCLUDE (ScreamUtils)

SET (num_mpi_ranks 1)
SET (config_defines "")

# Test utilities (f90)
CreateUnitTest(util_f90 util_test.f90 scream_share num_mpi_ranks config_defines)

# Test utilities (c++)
CreateUnitTest(util_cxx util_tests.cpp scream_share num_mpi_ranks config_defines)

# Test packs
CreateUnitTest(packs pack_tests.cpp scream_share num_mpi_ranks config_defines)
