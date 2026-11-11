# Homme system tests: theta-l DCMIP 2012 tests.
# See test/README.md for how tests are declared, and cmake/HommeTesting.cmake for the functions' docs.

homme_add_system_test (thetah-test22
  EXEC          theta-l-nlev20
  RANKS         16
  NAMELISTS     thetah-test22.nl
  OUTPUTS       dcmip2012_test2_21.nc)

homme_add_system_test (thetanh-test22
  EXEC          theta-l-nlev20
  RANKS         16
  NAMELISTS     thetanh-test22.nl
  OUTPUTS       dcmip2012_test2_21.nc)

# F90 and C++ versions, for each testing profile (see homme_set_profile_params)
if (HOMMEXX_BFB_TESTING)
  foreach (profile IN LISTS HOMME_TESTING_PROFILES)
    homme_set_profile_params (theta ${profile})
    set (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
    set (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
    homme_add_f90_cxx_system_tests (theta-fdc12-test21 CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC          theta-nlev128
      RANKS         16
      NAMELISTS     theta-fdc12-test21.nl
      VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/sab*-128.ascii
      OUTPUTS       dcmip2012_test2_11.nc)

    set (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
    set (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
    homme_add_f90_cxx_system_tests (theta-fdc12-test22 CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC          theta-nlev128
      RANKS         16
      NAMELISTS     theta-fdc12-test22.nl
      VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/sab*-128.ascii
      OUTPUTS       dcmip2012_test2_21.nc)

    set (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
    set (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
    homme_add_f90_cxx_system_tests (theta-fdc12-test3 CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC          theta-nlev128
      RANKS         16
      NAMELISTS     theta-fdc12-test3.nl
      VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/sab*-128.ascii
      OUTPUTS       dcmip2012_test31.nc)
  endforeach()
endif()
