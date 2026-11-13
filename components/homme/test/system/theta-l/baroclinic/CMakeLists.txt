# Homme system tests: theta-l baroclinic wave tests.
# See test/README.md for how tests are declared, and cmake/HommeTesting.cmake for the functions' docs.

homme_add_system_test (thetah-TC
  EXEC          theta-l-nlev30
  RANKS         16
  NAMELISTS     thetah-TC.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*30*
  OUTPUTS       dcmip2016_test21.nc)

homme_add_system_test (thetanh-TC
  EXEC          theta-l-nlev30
  RANKS         16
  NAMELISTS     thetanh-TC.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*30*
  OUTPUTS       dcmip2016_test21.nc)

homme_add_system_test (thetanh-TC-nudiv
  EXEC          theta-l-nlev30
  RANKS         16
  NAMELISTS     thetanh-TC-nudiv.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*30*
  OUTPUTS       dcmip2016_test21.nc)

homme_add_system_test (thetanh-TC-nutop
  EXEC          theta-l-nlev30
  RANKS         16
  NAMELISTS     thetanh-TC-nutop.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*30*
  OUTPUTS       dcmip2016_test21.nc)

homme_add_system_test (thetanhwet-TC
  EXEC          theta-l-nlev30
  RANKS         16
  NAMELISTS     thetanhwet-TC.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*30*
  OUTPUTS       dcmip2016_test21.nc)

homme_add_system_test (thetah-TC-ftype4
  EXEC          theta-l-nlev30
  RANKS         16
  NAMELISTS     thetah-TC-ftype4.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*30*
  OUTPUTS       dcmip2016_test21.nc)

# Jablonowski-Williamson baroclinic wave tests with various options, F90 and C++, for each
# testing profile (see homme_set_profile_params), used to check that F90 and C++ give the same answers.
# The test name describes the options: f=theta_advect_form, tt=ttype, hvs=hypervis_subcycle,
# hvst=hypervis_subcycle_tom, r=rsplit, qz=qsize, nutop=nu_top on/off, GB=GB tests (sl=semi-Lagrangian transport).
# NOTE: names are not parsed: the options are passed explicitly (and some deliberately differ from the name).
function (add_theta_baroclinic_tests name)
  cmake_parse_arguments (PARSE_ARGV 1 atbt "" "NAMELIST;MOISTURE;FORM;TTYPE;HY_MODE;HVS;HVST;RSPLIT;QSIZE;QSPLIT;NUTOP" "")
  set (HOMME_TEST_LIM 9)
  set (HOMME_TEST_MOISTURE ${atbt_MOISTURE})
  set (HOMME_TEST_HVSCALING 0) # const HV for now, tensor is tested in preqx
  set (HOMME_THETA_FORM ${atbt_FORM})
  set (HOMME_TTYPE ${atbt_TTYPE})
  set (HOMME_TEST_HVS ${atbt_HVS})
  set (HOMME_TEST_HVS_TOM ${atbt_HVST})
  set (HOMME_TEST_RSPLIT ${atbt_RSPLIT})
  set (HOMME_TEST_QSIZE ${atbt_QSIZE})
  set (HOMME_TEST_QSPLIT ${atbt_QSPLIT})
  set (HOMME_TEST_NUTOP ${atbt_NUTOP})
  set (HOMME_THETA_HY_MODE ${atbt_HY_MODE})
  set (HOMME_TEST_TIME_STEP 600)
  set (HOMME_TEST_VCOORD_INT_FILE sabi-128.ascii)
  set (HOMME_TEST_VCOORD_MID_FILE sabm-128.ascii)
  foreach (profile IN LISTS HOMME_TESTING_PROFILES)
    homme_set_profile_params (theta ${profile})
    homme_add_f90_cxx_system_tests (${name} CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC         theta-nlev128
      RANKS        16
      NAMELISTS    ${atbt_NAMELIST}
      VCOORD_FILES ${HOMME_TEST_DATA_DIR}/vcoord/sab*-128.ascii
      OUTPUTS      jw_baroclinic1.nc jw_baroclinic2.nc)
  endforeach()
endfunction()

if (HOMMEXX_BFB_TESTING)
  add_theta_baroclinic_tests (theta-f0-tt5-hvs1-hvst0-r3-qz1-nutopoff
    NAMELIST theta.nl MOISTURE dry FORM 0 TTYPE 5 HY_MODE true
    HVS 1 HVST 0 RSPLIT 3 QSIZE 1 QSPLIT "" NUTOP 0)
  add_theta_baroclinic_tests (theta-f1-tt5-hvs1-hvst0-r3-qz1-nutopoff
    NAMELIST theta.nl MOISTURE dry FORM 1 TTYPE 5 HY_MODE true
    HVS 1 HVST 0 RSPLIT 3 QSIZE 1 QSPLIT "" NUTOP 0)
  add_theta_baroclinic_tests (theta-f1-tt5-hvs1-hvst0-r0-qz1-nutopoff
    NAMELIST theta.nl MOISTURE dry FORM 1 TTYPE 5 HY_MODE true
    HVS 1 HVST 0 RSPLIT 0 QSIZE 1 QSPLIT "" NUTOP 0)
  add_theta_baroclinic_tests (theta-f1-tt10-hvs3-hvst0-r3-qz1-nutopoff
    NAMELIST theta.nl MOISTURE dry FORM 1 TTYPE 10 HY_MODE false
    HVS 3 HVST 0 RSPLIT 3 QSIZE 1 QSPLIT "" NUTOP 0)
  add_theta_baroclinic_tests (theta-f1-tt10-hvs3-hvst5-r3-qz1-nutopon
    NAMELIST theta.nl MOISTURE dry FORM 2 TTYPE 10 HY_MODE false
    HVS 3 HVST 5 RSPLIT 3 QSIZE 1 QSPLIT "" NUTOP 2.5e5)
  add_theta_baroclinic_tests (theta-f1-tt10-hvs1-hvst0-r2-qz10-nutopoff-GB
    NAMELIST theta.nl MOISTURE dry FORM 1 TTYPE 10 HY_MODE false
    HVS 1 HVST 0 RSPLIT 2 QSIZE 10 QSPLIT "" NUTOP 0)
  add_theta_baroclinic_tests (theta-f1-tt10-hvs1-hvst0-r2-qz10-nutopoff-GB-sl
    NAMELIST theta-sl.nl MOISTURE wet FORM 1 TTYPE 10 HY_MODE false
    HVS 1 HVST 0 RSPLIT 2 QSIZE 10 QSPLIT 6 NUTOP 0)
endif()
