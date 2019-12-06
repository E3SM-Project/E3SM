set (pfunit_sources
  test_map_lnd2rof_irrig.pf
  )

create_pFUnit_test(map_lnd2rof_irrig map_lnd2rof_irrig_exe
  "${pfunit_sources}" "")

target_link_libraries(map_lnd2rof_irrig_exe ${DRV_UNIT_TEST_LIBS})
