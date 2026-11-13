# Homme system tests: hydrostatic primitive equations (preqx and preqx_kokkos).
# See test/README.md for how tests are declared, and cmake/HommeTesting.cmake for the functions' docs.

homme_add_system_test (baro2b
  EXEC          baroC
  RANKS         16
  NAMELISTS     baro2b.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*26*
  OUTPUTS       baro2b-asp_baroclinic1.nc baro2b-asp_baroclinic2.nc)

homme_add_system_test (baro2c
  EXEC          baroC
  RANKS         16
  NAMELISTS     baro2c-run1.nl baro2c-run2.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*26*
  OUTPUTS       baro2c-run1-jw_baroclinic1.nc baro2c-run2-jw_baroclinic-000000000-1.nc
  OMP_NAMELISTS baro2c-run2-omp.nl
  OMP_THREADS   4
  REF_OUTPUTS   baro2c-run2-jw_baroclinic-000000000-1.nc
  CHECK_OUTPUTS baro2c-run2-omp-jw_baroclinic-000000000-1.nc)

homme_add_system_test (baro2d
  EXEC          baroC
  RANKS         16
  NAMELISTS     baro2d.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*26*
  OUTPUTS
    baro2d-asp_baroclinic1.nc
    baro2d-asp_baroclinic2.nc
    baro2d-asp_baroclinic3.nc)

homme_add_system_test (baroCamMoist
  EXEC          baroCam
  RANKS         16
  NAMELISTS     baroCamMoist.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*26*
  OUTPUTS       camBaroMoist-asp_baroclinic1.nc camBaroMoist-asp_baroclinic2.nc
  OMP_NAMELISTS baroCamMoist-omp3.nl baroCamMoist-omp4.nl
  OMP_THREADS   4
  REF_OUTPUTS
    camBaroMoist-asp_baroclinic1.nc
    camBaroMoist-asp_baroclinic2.nc
    camBaroMoist-asp_baroclinic1.nc
    camBaroMoist-asp_baroclinic2.nc
  CHECK_OUTPUTS
    camBaroMoist-omp3-asp_baroclinic1.nc
    camBaroMoist-omp3-asp_baroclinic2.nc
    camBaroMoist-omp4-asp_baroclinic1.nc
    camBaroMoist-omp4-asp_baroclinic2.nc)

homme_add_system_test (preqx-TC-ftype4
  EXEC          preqx-nlev30-interp
  RANKS         16
  NAMELISTS     preqx-TC-ftype4.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*30*
  OUTPUTS       dcmip2016_test21.nc)

# Nonhydrostatic gravity waves (planar), F90 and C++
homme_add_f90_cxx_system_tests (preqx-nhgw CXX_VS_F90
  EXEC          preqx-nlev10-native
  RANKS         4
  NAMELISTS     preqx-nhgw.nl
  OUTPUTS       planar_nonhydro_gravity_wave1.nc)

homme_add_f90_cxx_system_tests (preqx-nhgw-slice CXX_VS_F90
  EXEC          preqx-nlev10-native
  RANKS         4
  NAMELISTS     preqx-nhgw-slice.nl
  OUTPUTS       planar_nonhydro_gravity_wave1.nc)

# Jablonowski-Williamson baroclinic wave tests, F90 and C++, for each testing profile
# (see homme_set_profile_params), used to check that F90 and C++ give the same answers
if (HOMMEXX_BFB_TESTING)
  foreach (profile IN LISTS HOMME_TESTING_PROFILES)
    homme_set_profile_params (preqx ${profile})
    set (HOMME_TEST_LIM 8)
    set (HOMME_TEST_QSIZE 1)
    set (HOMME_TEST_RSPLIT 0)
    set (HOMME_TEST_MOISTURE dry)
    set (HOMME_TEST_HVSCALING 0)
    set (HOMME_TEST_NU 7e15)
    set (HOMME_TEST_NUDIV 7e15)
    set (HOMME_TEST_NUTOP 0)
    set (HOMME_TEST_TIME_STEP 600)
    set (HOMME_TEST_VCOORD_INT_FILE cami-26.ascii)
    set (HOMME_TEST_VCOORD_MID_FILE camm-26.ascii)
    homme_add_f90_cxx_system_tests (preqx-nlev26-dry-r0-samenu-consthv-lim8-q1 CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC          preqx-nlev26
      RANKS         16
      NAMELISTS     preqx.nl
      VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/cam*-26.ascii
      OUTPUTS       jw_baroclinic1.nc jw_baroclinic2.nc)

    homme_set_profile_params (preqx ${profile})
    set (HOMME_TEST_LIM 9)
    set (HOMME_TEST_QSIZE 6)
    set (HOMME_TEST_RSPLIT 3)
    set (HOMME_TEST_MOISTURE dry)
    set (HOMME_TEST_HVSCALING 0)
    set (HOMME_TEST_NU 7e15)
    set (HOMME_TEST_NUDIV 1e15)
    set (HOMME_TEST_NUTOP 2.5e5)
    set (HOMME_TEST_TIME_STEP 600)
    set (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
    set (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)
    homme_add_f90_cxx_system_tests (preqx-nlev72-dry-r3-diffnu-consthv-lim9-q6 CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC          preqx-nlev72
      RANKS         16
      NAMELISTS     preqx.nl
      VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/acme-72*.ascii
      OUTPUTS       jw_baroclinic1.nc jw_baroclinic2.nc)

    homme_set_profile_params (preqx-tensorhv ${profile})
    set (HOMME_TEST_LIM 9)
    set (HOMME_TEST_QSIZE 1)
    set (HOMME_TEST_RSPLIT 3)
    set (HOMME_TEST_MOISTURE moist)
    set (HOMME_TEST_HVSCALING 3.2)
    set (HOMME_TEST_NU 4.9e-7)
    set (HOMME_TEST_NUDIV 4.9e-7)
    set (HOMME_TEST_NUTOP 0)
    set (HOMME_TEST_TIME_STEP 300)
    set (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
    set (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)
    homme_add_f90_cxx_system_tests (preqx-nlev72-moist-r3-samenu-tensorhv-lim9-q1 CXX_VS_F90 NAME_SUFFIX ${HOMME_PROFILE_SUFFIX} LABELS ${profile}
      EXEC          preqx-nlev72
      RANKS         16
      NAMELISTS     preqx.nl
      VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/acme-72*.ascii
      OUTPUTS       jw_baroclinic1.nc jw_baroclinic2.nc)
  endforeach()
endif()
