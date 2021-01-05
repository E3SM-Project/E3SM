INCLUDE (ScreamUtils)

SET (NEED_LIBS shoc physics_share scream_share)
set(SHOC_TESTS_SRCS
    shoc_tests.cpp
    shoc_grid_tests.cpp
    shoc_diag_obklen_tests.cpp
    shoc_impli_comp_tmpi_tests.cpp
    shoc_impli_dp_inverse_tests.cpp
    shoc_impli_sfc_fluxes_tests.cpp
    shoc_impli_srf_stress_tests.cpp
    shoc_impli_srf_tke_tests.cpp
    shoc_vertflux_tests.cpp
    shoc_varorcovar_tests.cpp
    shoc_energy_fixer_tests.cpp
    shoc_energy_update_dse_tests.cpp
    shoc_energy_integral_tests.cpp
    shoc_energy_total_fixer_tests.cpp
    shoc_energy_dse_fixer_tests.cpp
    shoc_energy_threshold_fixer_tests.cpp
    shoc_length_tests.cpp
    shoc_brunt_length_tests.cpp
    shoc_l_inf_length_tests.cpp
    shoc_conv_vel_length_tests.cpp
    shoc_check_length_tests.cpp
    shoc_conv_time_length_tests.cpp
    shoc_mix_length_tests.cpp
    shoc_fterm_input_third_moms_tests.cpp
    shoc_fterm_diag_third_moms_tests.cpp
    shoc_omega_diag_third_moms_tests.cpp
    shoc_xy_diag_third_moms_tests.cpp
    shoc_aa_diag_third_moms_tests.cpp
    shoc_w3_diag_third_moms_tests.cpp
    shoc_clip_third_moms_tests.cpp
    shoc_tke_tests.cpp
    shoc_check_tke_tests.cpp
    shoc_tke_column_stab_tests.cpp
    shoc_tke_shr_prod_tests.cpp
    shoc_tke_isotropic_ts_tests.cpp
    shoc_tke_adv_sgs_tke_tests.cpp
    shoc_eddy_diffusivities_tests.cpp
    shoc_diag_second_mom_srf_test.cpp
    shoc_linear_interp_tests.cpp
    shoc_diag_third_tests.cpp
    shoc_compute_diag_third_tests.cpp
    shoc_assumed_pdf_tests.cpp
    shoc_pdf_tildetoreal_tests.cpp
    shoc_pdf_vv_parameters_tests.cpp
    shoc_pdf_thl_parameters_tests.cpp
    shoc_pdf_qw_parameters_tests.cpp
    shoc_pdf_inplume_corr_tests.cpp
    shoc_pdf_computetemp_tests.cpp
    shoc_pdf_compute_qs_tests.cpp
    shoc_pdf_compute_s_tests.cpp
    shoc_pdf_compute_buoyflux_tests.cpp
    shoc_pdf_compute_liqflux_tests.cpp
    shoc_pdf_compute_cloudvar_tests.cpp
    shoc_pdf_compute_sgsliq_tests.cpp
    shoc_diag_second_mom_ubycond_test.cpp
    shoc_pblintd_init_pot_test.cpp
    shoc_unit_tests.cpp
    shoc_diag_second_moments_lbycond_tests.cpp
    shoc_diag_second_moments_tests.cpp
    shoc_diag_second_shoc_moments_tests.cpp
    shoc_pblintd_cldcheck_tests.cpp
    shoc_compute_shoc_vapor_tests.cpp
    shoc_update_prognostics_implicit_tests.cpp
    shoc_main_tests.cpp
    shoc_pblintd_height_tests.cpp
    shoc_vd_shoc_decomp_and_solve_tests.cpp
    shoc_pblintd_surf_temp_tests.cpp
    shoc_pblintd_check_pblh_tests.cpp
    shoc_pblintd_tests.cpp
    ) # SHOC_TESTS_SRCS

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  CreateUnitTest(shoc_tests "${SHOC_TESTS_SRCS}" "${NEED_LIBS}" THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC} DEP shoc_tests_ut_np1_omp1)
endif()

CreateUnitTest(shoc_run_and_cmp "shoc_run_and_cmp.cpp" "${NEED_LIBS}" THREADS ${SCREAM_TEST_MAX_THREADS} EXE_ARGS "${SCREAM_TEST_DATA_DIR}/shoc_run_and_cmp.baseline" DEP shoc_tests_ut_np1_omp1 EXCLUDE_MAIN_CPP)

# By default, baselines should be created using all fortran
add_custom_target(shoc_baseline
  COMMAND $<TARGET_FILE:shoc_run_and_cmp> -f -g ${SCREAM_TEST_DATA_DIR}/shoc_run_and_cmp.baseline)

add_dependencies(baseline shoc_baseline)
