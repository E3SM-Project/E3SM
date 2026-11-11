# Homme system tests: theta-l Held-Suarez tests.
# See test/README.md for how tests are declared, and cmake/HommeTesting.cmake for the functions' docs.

# F90 and C++ versions, for each testing profile (see homme_set_profile_params)
#   theta-fhs1: no topo, moist, ftype0, EUL, tensor HV, NH, ttype10, v_alg1
#   theta-fhs2: no topo, moist, ftype0, SL,  tensor HV, NH, ttype10, v_alg10
#   theta-fhs3: no topo, moist, ftype2, SL,  tensor HV, NH, ttype9
if (HOMMEXX_BFB_TESTING)
  foreach (profile IN LISTS HOMME_TESTING_PROFILES)
    homme_set_profile_params (theta ${profile})
    set (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
    set (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
    homme_add_f90_cxx_system_tests (theta-fhs1 CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC          theta-nlev128
      RANKS         16
      NAMELISTS     thetahs1.nl
      VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/sab*-128.ascii
      OUTPUTS       held_suarez01.nc held_suarez02.nc)

    set (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
    set (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
    homme_add_f90_cxx_system_tests (theta-fhs2 CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC          theta-nlev128
      RANKS         16
      NAMELISTS     thetahs2.nl
      VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/sab*-128.ascii
      OUTPUTS       held_suarez01.nc held_suarez02.nc)

    set (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
    set (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
    homme_add_f90_cxx_system_tests (theta-fhs3 CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC          theta-nlev128
      RANKS         16
      NAMELISTS     thetahs3.nl
      VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/sab*-128.ascii
      OUTPUTS       held_suarez01.nc held_suarez02.nc)
  endforeach()
endif()
