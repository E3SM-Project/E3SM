set(pfunit_sources
  test_update_landunit_weights_one_gcell.pf
  test_update_landunit_weights.pf)

create_pFUnit_test(dynLandunitArea test_dynLandunitArea_exe
  "${pfunit_sources}" "")

target_link_libraries(test_dynLandunitArea_exe clm csm_share)

