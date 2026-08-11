include(ScreamUtils)

add_subdirectory(infra)

set(ZM_TESTS_SRCS
  zm_ientropy_tests.cpp
  zm_entropy_tests.cpp
  zm_transport_tracer_tests.cpp
  zm_transport_momentum_tests.cpp
  zm_compute_dilute_cape_tests.cpp
  zm_find_mse_max_tests.cpp
  zm_compute_dilute_parcel_tests.cpp
  zm_compute_cape_from_parcel_tests.cpp
  zm_conv_mcsp_calculate_shear_tests.cpp
  zm_conv_mcsp_tend_tests.cpp
  zm_conv_main_tests.cpp
  zm_conv_evap_tests.cpp
  zm_calc_fractional_entrainment_tests.cpp
  zm_downdraft_properties_tests.cpp
  zm_cloud_properties_tests.cpp
  zm_closure_tests.cpp
  zm_calc_output_tend_tests.cpp
) # ZM_TESTS_SRCS

CreateUnitTest(zm_tests
  SOURCES ${ZM_TESTS_SRCS}
  LIBS zm zm_test_infra
  EXE_ARGS "--args ${SCREAM_BASELINE_FILE_ARG}"
  THREADS ${SCREAM_BASELINE_TEST_THREADS}
  LABELS zm physics baseline_gen baseline_cmp
)
