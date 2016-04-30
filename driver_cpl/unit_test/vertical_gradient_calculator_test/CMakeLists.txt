set (pfunit_sources
  test_vertical_gradient_calculator_2nd_order.pf
  test_vertical_gradient_calculator_continuous.pf
  test_vertical_gradient_calculator_factory.pf
  )

# Only include output_vertical_gradient_calculator_continuous.pf if the
# environment variable OUTPUT_VGCONT is set: The "tests" in here aren't true
# unit tests, but rather more like "functional unit tests" (printing output for
# various inputs, for later inspection). This environment variable mechanism is
# a kludge to prevent running these "tests" that do i/o most of the time, only
# including them when requested.
#
# However, this mechanism requires rerunning cmake if you set or unset this
# environment variable.
if (DEFINED ENV{OUTPUT_VGCONT})
  list(APPEND pfunit_sources
    output_vertical_gradient_calculator_continuous.pf
    )
endif()

set (extra_sources
  vertical_gradient_calculator_continuousNoLimit.F90
  )

create_pFUnit_test(vertical_gradient_calculator vertical_gradient_calculator_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(vertical_gradient_calculator_exe ${DRV_UNIT_TEST_LIBS})
