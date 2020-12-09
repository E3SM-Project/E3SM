# Lists of test files for the HOMME regression tests
#### WIP

#proposed tests:
#nlev26-dry-rsplit0-(nudiv==nu)-consthv-lim8-qsplit1
#nlev72-moist-rsplit3-(nudiv!=nu)-consthv-lim8-qsplit6
#nlev72-moist-rsplit3-(nudiv==nu)-tensorhv-lim9-qsplit6

#//this is what is tested now -- 26 and 72 levels,  dry/moist, rsplit 0 and 3, 
#nudiv==nu and nudiv != nu, tensorhv and const hv, lim9 and lim8 (do we want lim8 still?), qsplit 1 and 6

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
  baroCamMoist-acc.cmake
  baro2d-imp.cmake
  thetah-test22.cmake
  thetanh-test22.cmake
  preqx-TC-ftype4.cmake
  thetah-TC-ftype4.cmake
  thetah-TC.cmake
  thetanh-TC.cmake
  thetanh-TC-nudiv.cmake
  thetanh-TC-nutop.cmake
  thetanhwet-TC.cmake
  hommetool.cmake
  thetah-nhgw.cmake
  thetanh-nhgw.cmake
  preqx-nhgw.cmake
  thetah-nhgw-slice.cmake
  thetanh-nhgw-slice.cmake
  preqx-nhgw-slice.cmake
  sweqx-dbvor.cmake
)

IF (${HOMME_ENABLE_COMPOSE})
  LIST(APPEND HOMME_TESTS
    thetah-sl-test11conv-r1t2-cdr20.cmake
    thetah-sl-test11conv-r0t1-cdr30-rrm.cmake
    thetah-sl-dcmip16_test1pg2.cmake
    )
ENDIF()

IF (${BUILD_HOMME_PREQX_KOKKOS})
  # Lists of test files for the HOMME kokkos regression tests.
  # These tests come in pairs, so that, besides checking
  # performance, we can compare the kokkos output with the
  # fortran one later on.
  SET(HOMME_PREQX_TESTS_WITH_PROFILE
    preqx-nlev26-dry-r0-samenu-consthv-lim8-q1.cmake
    preqx-nlev26-dry-r0-samenu-consthv-lim8-q1-kokkos.cmake
    preqx-nlev72-dry-r3-diffnu-consthv-lim9-q6.cmake
    preqx-nlev72-dry-r3-diffnu-consthv-lim9-q6-kokkos.cmake
    preqx-nlev72-moist-r3-samenu-tensorhv-lim9-q1.cmake
    preqx-nlev72-moist-r3-samenu-tensorhv-lim9-q1-kokkos.cmake
  )

  #This list (COMPARE_F_C_TEST) contains tests for which
  #F vc C comparison will be run.
  #Note: we run both a cprnc on the output nc files AND
  #      a comparison of the values of diagnostic quantities
  #      on the raw output files
  IF (${ENABLE_PREQX_KOKKOS_BFB_TESTS})
    SET (PREQX_COMPARE_F_C_TEST
    preqx-nlev26-dry-r0-samenu-consthv-lim8-q1
    preqx-nlev72-dry-r3-diffnu-consthv-lim9-q6
    preqx-nlev72-moist-r3-samenu-tensorhv-lim9-q1
    )
  ENDIF ()
ENDIF()
