# Lists of test files for the HOMME regression tests
SET(HOMME_TESTS 
  swtc1.cmake
  swtc2.cmake
  swtc5.cmake
  swtc6.cmake
  swimtc5.cmake
  baro2b.cmake
  baro2c.cmake
  baro2d.cmake
  baroCamMoist.cmake
  baroCamMoist-SL.cmake
  baroCamMoist-acc.cmake
  baro2d-imp.cmake
  thetah-test22.cmake
  thetanh-test22.cmake
  thetah-TC.cmake
  thetanh-TC.cmake
  thetanhwet-TC.cmake
  templates.cmake
)

IF (${BUILD_HOMME_PREQX_KOKKOS})
  # Lists of test files for the HOMME kokkos regression tests.
  # These tests come in pairs, so that, besides checking
  # performance, we can compare the kokkos output with the
  # fortran one later on.
  SET(HOMME_PREQX_TESTS_WITH_PROFILE
    preqx-nlev26-qsize4-r0-dry.cmake
    preqx-nlev26-qsize4-r0-dry-kokkos.cmake
    preqx-nlev26-qsize4-r3-dry.cmake
    preqx-nlev26-qsize4-r3-dry-kokkos.cmake
    preqx-nlev26-qsize4-r0-moist.cmake
    preqx-nlev26-qsize4-r0-moist-kokkos.cmake
    preqx-nlev26-qsize4-r3-moist.cmake
    preqx-nlev26-qsize4-r3-moist-kokkos.cmake
    preqx-nlev72-qsize4-r0-dry.cmake
    preqx-nlev72-qsize4-r0-dry-kokkos.cmake
    preqx-nlev72-qsize4-r3-dry.cmake
    preqx-nlev72-qsize4-r3-dry-kokkos.cmake
    preqx-nlev72-qsize4-r0-moist.cmake
    preqx-nlev72-qsize4-r0-moist-kokkos.cmake
    preqx-nlev72-qsize4-r3-moist.cmake
    preqx-nlev72-qsize4-r3-moist-kokkos.cmake
    preqx-nlev72-qsize4-r3-q6-dry.cmake
    preqx-nlev72-qsize4-r3-q6-dry-kokkos.cmake
    preqx-nlev72-qsize4-r3-tensorhv-dry.cmake
    preqx-nlev72-qsize4-r3-tensorhv-dry-kokkos.cmake
    preqx-nlev72-qsize4-r3-nudiv-dry.cmake
    preqx-nlev72-qsize4-r3-nudiv-dry-kokkos.cmake
    preqx-nlev72-qsize10-r3-lim9-dry.cmake
    preqx-nlev72-qsize10-r3-lim9-dry-kokkos.cmake
  )

  #This list (COMPARE_F_C_TEST) contains tests for which
  #F vc C comparison will be run.
  #Note: we run both a cprnc on the output nc files AND
  #      a comparison of the values of diagnostic quantities
  #      on the raw output files
  SET (PREQX_COMPARE_F_C_TEST
    preqx-nlev26-qsize4-r0-dry
    preqx-nlev26-qsize4-r3-dry
    preqx-nlev26-qsize4-r0-moist
    preqx-nlev26-qsize4-r3-moist
    preqx-nlev72-qsize4-r0-dry
    preqx-nlev72-qsize4-r3-dry
    preqx-nlev72-qsize4-r0-moist
    preqx-nlev72-qsize4-r3-moist
    preqx-nlev72-qsize4-r3-q6-dry
    preqx-nlev72-qsize4-r3-tensorhv-dry
    preqx-nlev72-qsize4-r3-nudiv-dry
    preqx-nlev72-qsize10-r3-lim9-dry
  )
ENDIF()
