set (pfunit_sources
  test_vertical_gradient_test_utils.pf
  test_vertical_gradient_calculator_2nd_order.pf
  test_vertical_gradient_calculator_continuous.pf
  )

set (extra_sources
  vertical_gradient_test_utils.F90
  )

create_pFUnit_test(vertical_gradient_calculator vertical_gradient_calculator_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(vertical_gradient_calculator_exe ${DRV_UNIT_TEST_LIBS})
