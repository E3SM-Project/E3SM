include (ScreamUtils)

add_subdirectory(infra)

set(SHOC_TESTS_SRCS
    shoc_tests.cpp
    shoc_grid_tests.cpp
    shoc_diag_obklen_tests.cpp
    shoc_impli_comp_tmpi_tests.cpp
    shoc_impli_dp_inverse_tests.cpp
    shoc_vertflux_tests.cpp
    shoc_varorcovar_tests.cpp
    shoc_energy_fixer_tests.cpp
    shoc_energy_update_dse_tests.cpp
    shoc_energy_integral_tests.cpp
    shoc_length_tests.cpp
    shoc_brunt_length_tests.cpp
    shoc_l_inf_length_tests.cpp
    shoc_check_length_tests.cpp
    shoc_mix_length_tests.cpp
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
    shoc_compute_shoc_temperature_tests.cpp
    ) # SHOC_TESTS_SRCS

# All tests should understand the same baseline args
if (SCREAM_ENABLE_BASELINE_TESTS)
  if (SCREAM_ONLY_GENERATE_BASELINES)
    set(BASELINE_FILE_ARG "-g -b ${SCREAM_BASELINES_DIR}/data")
    # We don't want to do thread spreads when generating. That
    # could cause race conditions in the file system.
    set(SHOC_THREADS "${SCREAM_TEST_MAX_THREADS}")
  else()
    set(BASELINE_FILE_ARG "-c -b ${SCREAM_BASELINES_DIR}/data")
    set(SHOC_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
  endif()
else()
  set(BASELINE_FILE_ARG "-n") # no baselines
  set(SHOC_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
endif()

CreateUnitTest(shoc_tests "${SHOC_TESTS_SRCS}"
  LIBS shoc shoc_test_infra
  THREADS ${SHOC_THREADS}
  EXE_ARGS "--args ${BASELINE_FILE_ARG}"
  LABELS "shoc;physics;baseline_gen;baseline_cmp"
  )

# If small kernels are ON, we don't need a separate executable to test them.
# Also, we never want to generate baselines with this separate executable
if (NOT SCREAM_SHOC_SMALL_KERNELS AND NOT SCREAM_ONLY_GENERATE_BASELINES)
  CreateUnitTest(shoc_sk_tests "shoc_main_tests.cpp"
    LIBS shoc_sk shoc_test_infra
    THREADS ${SHOC_THREADS}
    EXE_ARGS "--args ${BASELINE_FILE_ARG}"
    LABELS "shoc;physics;baseline_cmp"
    )
endif()

CreateUnitTest(shoc_run_and_cmp "shoc_run_and_cmp.cpp"
  LIBS shoc shoc_test_infra
  EXCLUDE_MAIN_CPP
  THREADS ${SCREAM_TEST_MAX_THREADS}
  EXE_ARGS "${BASELINE_FILE_ARG}"
  LABELS "shoc;physics;baseline_gen;baseline_cmp")
