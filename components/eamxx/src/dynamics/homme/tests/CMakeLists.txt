# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT SCREAM_BASELINES_ONLY)
  include (ScreamUtils)

  # Get or create the dynamics lib
  #                 HOMME_TARGET   NP PLEV QSIZE_D
  CreateDynamicsLib("theta-l_kokkos"  4   72   10)

  # Test dynamics-physics fields remapping
  CreateUnitTest(homme_pd_remap
      "homme_pd_remap_tests.cpp;test_helper_mod.F90"
      "scream_share;${dynLibName}"
      MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
      LABELS "dynamics")

  # Test I/O on dyn grid
  CreateUnitTest(dyn_grid_io
      "dyn_grid_io.cpp;test_helper_mod.F90"
      "scream_share;scream_io;${dynLibName}"
      MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
      LABELS "dynamics;io")
endif()
