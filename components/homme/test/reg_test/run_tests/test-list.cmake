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
  thetanh-moist-bubble.cmake
  thetanh-dry-bubble.cmake
)

IF (HOMME_ENABLE_COMPOSE)
  LIST(APPEND HOMME_TESTS
    thetah-sl-test11conv-r1t2-cdr20.cmake
    thetah-sl-test11conv-r0t1-cdr30-rrm.cmake
    thetah-sl-dcmip16_test1pg2.cmake
    thetah-sl-testconv-3e.cmake)
ENDIF()

SET(HOMME_RUN_TESTS_DIR ${HOMME_SOURCE_DIR}/test/reg_test/run_tests)

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
  IF (${HOMMEXX_BFB_TESTING})
    SET (PREQX_COMPARE_F_C_TEST
    preqx-nlev26-dry-r0-samenu-consthv-lim8-q1
    preqx-nlev72-dry-r3-diffnu-consthv-lim9-q6
    preqx-nlev72-moist-r3-samenu-tensorhv-lim9-q1
    )
  ENDIF ()

  LIST(APPEND HOMME_TESTS
    preqx-nhgw-kokkos.cmake
    preqx-nhgw-slice-kokkos.cmake)
  IF (HOMMEXX_BFB_TESTING)
    LIST(APPEND HOMME_ONEOFF_CVF_TESTS
      preqx-nhgw
      preqx-nhgw-slice)
  ENDIF()
ENDIF()

IF (BUILD_HOMME_THETA_KOKKOS)
  # Various one-off tests.
  IF (HOMME_ENABLE_COMPOSE)
    LIST(APPEND HOMME_TESTS
      thetah-sl-test11conv-r0t1-cdr30-rrm-kokkos.cmake
      thetah-sl-testconv-3e-kokkos.cmake)
    IF (HOMMEXX_BFB_TESTING)
      LIST(APPEND HOMME_ONEOFF_CVF_TESTS
        thetah-sl-test11conv-r0t1-cdr30-rrm)
    ENDIF()
  ENDIF()
  LIST(APPEND HOMME_TESTS
    thetanh-moist-bubble-kokkos.cmake
    thetanh-dry-bubble-kokkos.cmake
    thetah-nhgw-kokkos.cmake
    thetanh-nhgw-kokkos.cmake
    thetah-nhgw-slice-kokkos.cmake
    thetanh-nhgw-slice-kokkos.cmake
    thetanh-moist-bubble-sl.cmake
    thetanh-moist-bubble-sl-kokkos.cmake
    thetanh-moist-bubble-sl-pg2.cmake
    thetanh-moist-bubble-sl-pg2-kokkos.cmake)
  IF (HOMMEXX_BFB_TESTING)
    LIST(APPEND HOMME_ONEOFF_CVF_TESTS
      thetanh-moist-bubble
      thetanh-dry-bubble
      thetah-nhgw
      thetanh-nhgw
      thetah-nhgw-slice
      thetanh-nhgw-slice
      thetanh-moist-bubble-sl
      thetanh-moist-bubble-sl-pg2)
  ENDIF()

  #cmake/namelist will be built with create-... script
  SET(TESTS_WITH_AUTO
     theta-f0-tt5-hvs1-hvst0-r3-qz1-nutopoff
     theta-f1-tt5-hvs1-hvst0-r3-qz1-nutopoff
     theta-f1-tt5-hvs1-hvst0-r0-qz1-nutopoff
     theta-f1-tt10-hvs3-hvst0-r3-qz1-nutopoff
     theta-f1-tt10-hvs3-hvst5-r3-qz1-nutopon
     theta-f1-tt10-hvs1-hvst0-r2-qz10-nutopoff-GB
  )
  
  #manually added, create namelists and cmake files for F and cxx by hand
  #all have to be named theta-f* to make setup work (this could be fixed)
#     theta-fhs1 Held-Suarez, no topo, moist, ftype0, EUL, tensor HV, NH, ttype10, v_alg1
#     theta-fhs2              no topo, moist, ftype0, SL,  tensor HV, NH, ttype10, v_agl10
#     theta-fhs3              no topo, moist, ftype2, SL,  tensor HV, NH,  ttype9
#     theta-fdc12-test21      ttype10
#     theta-fdc12-test22      ttype9
#     theta-fdc12-test3       ttype9


  SET(TESTS_WITHOUT_AUTO
     theta-f1-tt10-hvs1-hvst0-r2-qz10-nutopoff-GB-sl
     theta-fhs1
     theta-fhs2
     theta-fhs3
     theta-fdc12-test21
     theta-fdc12-test22
     theta-fdc12-test3
  )

  #all tests that will be used for cxx-vs-F bfb testing
  SET(HOMME_THETA_TESTS_WITH_PROFILE_1
     ${TESTS_WITHOUT_AUTO} ${TESTS_WITH_AUTO}
  )

  set(HOMME_THETA_TESTS_WITH_PROFILE "")
  FOREACH(JJ ${HOMME_THETA_TESTS_WITH_PROFILE_1})
    LIST(APPEND HOMME_THETA_TESTS_WITH_PROFILE
        ${JJ}.cmake
        ${JJ}-kokkos.cmake
    )
  ENDFOREACH()
  IF (HOMMEXX_BFB_TESTING)
    SET (THETA_COMPARE_F_C_TEST ${HOMME_THETA_TESTS_WITH_PROFILE_1})
  ENDIF ()
ENDIF()

IF (BUILD_HOMME_PREQX_KOKKOS OR BUILD_HOMME_THETA_KOKKOS)
  IF (HOMMEXX_BFB_TESTING)
    CREATE_CXX_VS_F90_TESTS_WITH_PROFILE(HOMME_ONEOFF_CVF_TESTS short)
  ENDIF()
ENDIF()
