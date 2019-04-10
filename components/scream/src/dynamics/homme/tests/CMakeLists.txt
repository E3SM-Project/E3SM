INCLUDE (ScreamUtils)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D USE_PIO USE_ENERGY
CreateDynamicsLib("preqx_kokkos"  4   72   4  FALSE   FALSE)

SET (EXTRA_INCLUDE_DIRS ${HOMME_INCLUDE_DIRS})
SET (EXTRA_CONFIG_DEFS
     HAVE_CONFIG_H
     HOMMEXX_CONFIG_IS_CMAKE
)

SET (NEED_LIBS ${dynLibName} scream_share pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
SET (UNIT_TEST_EXTRA_LIBRARY_DIRS ${dynLibDir})

# Test dynamics-physics fields remapping
CreateUnitTest(homme_pd_remap "homme_pd_remap_tests.cpp" "${NEED_LIBS}" CONFIG_DEFS ${EXTRA_CONFIG_DEFS} INCLUDE_DIRS ${EXTRA_INCLUDE_DIRS})
