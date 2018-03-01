create_pFUnit_test(daylength test_daylength_exe
  "test_daylength.pf" "")

target_link_libraries(test_daylength_exe clm csm_share)