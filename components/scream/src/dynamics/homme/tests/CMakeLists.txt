# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT ${SCREAM_BASELINES_ONLY})
  include (ScreamUtils)

  # Get or create the dynamics lib
  #                 HOMME_TARGET   NP PLEV QSIZE_D USE_PIO USE_ENERGY
  CreateDynamicsLib("preqx_kokkos"  4   72   4  FALSE   FALSE)

  set (EXTRA_CONFIG_DEFS
       HAVE_CONFIG_H
       HOMMEXX_CONFIG_IS_CMAKE
  )

  set (NEED_LIBS ${dynLibName} scream_share pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  set (UNIT_TEST_EXTRA_LIBRARY_DIRS ${dynLibDir})

  # Test dynamics-physics fields remapping
  CreateUnitTest(homme_pd_remap "homme_pd_remap_tests.cpp" "${NEED_LIBS}"
    COMPILER_DEFS ${EXTRA_CONFIG_DEFS})
endif()
