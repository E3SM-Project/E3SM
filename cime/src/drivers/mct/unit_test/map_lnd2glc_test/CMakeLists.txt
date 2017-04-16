set (pfunit_sources
  test_map_lnd2glc.pf
  )

create_pFUnit_test(map_lnd2glc map_lnd2glc_exe
  "${pfunit_sources}" "")

target_link_libraries(map_lnd2glc_exe ${DRV_UNIT_TEST_LIBS})
