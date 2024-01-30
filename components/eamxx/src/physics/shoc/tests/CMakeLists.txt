include (ScreamUtils)

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
    shoc_check_length_tests.cpp
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
    shoc_compute_shoc_temperature_tests.cpp
    ) # SHOC_TESTS_SRCS

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  CreateUnitTest(shoc_tests    "${SHOC_TESTS_SRCS}"
    LIBS shoc
    THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
  )

  if (NOT SCREAM_SMALL_KERNELS)
    CreateUnitTest(shoc_sk_tests "${SHOC_TESTS_SRCS}"
      LIBS shoc_sk
      THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
      EXE_ARGS shoc_main_bfb
    )
  endif()
endif()

if (SCREAM_ENABLE_BASELINE_TESTS)
  set(BASELINE_FILE_ARG "-b ${SCREAM_TEST_DATA_DIR}/shoc_run_and_cmp.baseline")

  CreateUnitTestExec(shoc_run_and_cmp "shoc_run_and_cmp.cpp"
    LIBS shoc
    EXCLUDE_MAIN_CPP)

  CreateUnitTestFromExec(shoc_run_and_cmp_cxx shoc_run_and_cmp
    THREADS ${SCREAM_TEST_MAX_THREADS}
    EXE_ARGS "${BASELINE_FILE_ARG}"
    LABELS "shoc;physics")

  CreateUnitTestFromExec(shoc_run_and_cmp_f90 shoc_run_and_cmp
    THREADS ${SCREAM_TEST_MAX_THREADS}
    EXE_ARGS "-f ${BASELINE_FILE_ARG}"
    LABELS "shoc;physics")

  #
  # Use fake tests to generate shell commands to generate baselines
  #
  CreateUnitTestFromExec(shoc_baseline_f90_fake shoc_run_and_cmp
    THREADS ${SCREAM_TEST_MAX_THREADS}
    EXE_ARGS "-f -g ${BASELINE_FILE_ARG}"
    PROPERTIES DISABLED True)

  CreateUnitTestFromExec(shoc_baseline_cxx_fake shoc_run_and_cmp
    THREADS ${SCREAM_TEST_MAX_THREADS}
    EXE_ARGS "-g ${BASELINE_FILE_ARG}"
    PROPERTIES DISABLED True)

  if (SCREAM_TEST_MAX_THREADS GREATER 1)
    get_test_property(shoc_baseline_f90_fake_omp${SCREAM_TEST_MAX_THREADS} FULL_TEST_COMMAND SHOC_F90_GEN)
    get_test_property(shoc_baseline_cxx_fake_omp${SCREAM_TEST_MAX_THREADS} FULL_TEST_COMMAND SHOC_CXX_GEN)
  else()
    get_test_property(shoc_baseline_f90_fake FULL_TEST_COMMAND SHOC_F90_GEN)
    get_test_property(shoc_baseline_cxx_fake FULL_TEST_COMMAND SHOC_CXX_GEN)
  endif()

  if (SHOC_F90_GEN STREQUAL "NOTFOUND")
    message(FATAL_ERROR "Could not get FULL_TEST_COMMAND for shoc_baseline fake test")
  endif()

  separate_arguments(SHOC_F90_GEN_ARGS UNIX_COMMAND "${SHOC_F90_GEN}")
  separate_arguments(SHOC_CXX_GEN_ARGS UNIX_COMMAND "${SHOC_CXX_GEN}")

  add_custom_target(shoc_baseline_f90
    COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=${SCREAM_TEST_MAX_THREADS} ${SHOC_F90_GEN_ARGS})

  add_custom_target(shoc_baseline_cxx
    COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=${SCREAM_TEST_MAX_THREADS} ${SHOC_CXX_GEN_ARGS})

  # By default, baselines should be created using all fortran (make baseline). If the user wants
  # to use CXX to generate their baselines, they should use "make baseline_cxx".
  add_dependencies(baseline     shoc_baseline_f90)
  add_dependencies(baseline_cxx shoc_baseline_cxx)
endif()
