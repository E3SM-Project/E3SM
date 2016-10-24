set (pfunit_sources
  test_vertical_gradient_calculator_2nd_order.pf
  )

create_pFUnit_test(vertical_gradient_calculator vertical_gradient_calculator_exe
  "${pfunit_sources}" "")

target_link_libraries(vertical_gradient_calculator_exe ${DRV_UNIT_TEST_LIBS})
