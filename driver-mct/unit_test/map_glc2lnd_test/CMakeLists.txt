set (pfunit_sources
  test_map_glc2lnd.pf
  )

create_pFUnit_test(map_glc2lnd map_glc2lnd_exe
  "${pfunit_sources}" "")

target_link_libraries(map_glc2lnd_exe ${DRV_UNIT_TEST_LIBS})
