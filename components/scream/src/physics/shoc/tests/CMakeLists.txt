INCLUDE (ScreamUtils)

SET (NEED_LIBS shoc physics_share scream_control scream_share)
set(SHOC_TESTS_SRCS
    shoc_tests.cpp
    shoc_grid_tests.cpp
    shoc_vertflux_tests.cpp
    shoc_varorcovar_tests.cpp
    shoc_energy_update_dse_tests.cpp
    shoc_energy_integral_tests.cpp
    shoc_energy_total_fixer_tests.cpp
    shoc_energy_dse_fixer_tests.cpp
    shoc_energy_threshold_fixer_tests.cpp
    shoc_brunt_length_tests.cpp
    shoc_l_inf_length_tests.cpp
    shoc_conv_vel_length_tests.cpp
    shoc_check_length_tests.cpp
    shoc_conv_time_length_tests.cpp
    shoc_mix_length_tests.cpp
    shoc_tke_column_stab_tests.cpp
    shoc_tke_shr_prod_tests.cpp
    shoc_tke_isotropic_ts_tests.cpp
    shoc_tke_adv_sgs_tke_tests.cpp
    shoc_eddy_diffusivities_tests.cpp
    shoc_diag_second_mom_srf_test.cpp
    shoc_linearinterp_tests.cpp
    shoc_unit_tests.cpp)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  CreateUnitTest(shoc_tests "${SHOC_TESTS_SRCS}" "${NEED_LIBS}" THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC} DEP shoc_tests_ut_np1_omp1)
  # CreateUnitTest(shoc_stand_alone "shoc_stand_alone.cpp" "${NEED_LIBS}")
endif()

CreateUnitTest(shoc_run_and_cmp "shoc_run_and_cmp.cpp" "${NEED_LIBS}" THREADS ${SCREAM_TEST_MAX_THREADS} EXE_ARGS "${SCREAM_TEST_DATA_DIR}/shoc_run_and_cmp.baseline" DEP shoc_tests_ut_np1_omp1 EXCLUDE_MAIN_CPP)

# By default, baselines should be created using all fortran
add_custom_target(shoc_baseline
  COMMAND $<TARGET_FILE:shoc_run_and_cmp> -f -g ${SCREAM_TEST_DATA_DIR}/shoc_run_and_cmp.baseline)

add_dependencies(baseline shoc_baseline)

