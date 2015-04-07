set (pfunit_sources
  test_irrigation_deficit.pf
  test_irrigation_singlepatch.pf
  test_irrigation_multipatch.pf)

# extra sources used for this test, which are not .pf files
set (extra_sources
  IrrigationWrapperMod.F90)

create_pFUnit_test(irrigation test_irrigation_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_irrigation_exe clm csm_share)