# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include (ScreamUtils)

  # Get or create the dynamics lib
  #                 HOMME_TARGET   NP PLEV QSIZE_D
  CreateDynamicsLib("theta-l_kokkos"  4   72   10)

  # Test dynamics-physics fields remapping
  CreateUnitTest(homme_pd_remap
    SOURCES homme_pd_remap_tests.cpp test_helper_mod.F90
    LIBS ${dynLibName}
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    LABELS dynamics
  )

  # Compare (and time) the BoundaryExchange-based and the p2p-based phys->dyn remaps
  CreateUnitTest(homme_unique_grid_p2d
    SOURCES homme_unique_grid_p2d_tests.cpp test_helper_mod.F90
    LIBS ${dynLibName}
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    LABELS dynamics
  )
endif()
