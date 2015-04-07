create_pFUnit_test(humanstress test_humanstress_exe
  "test_humanstress.pf" "")

target_link_libraries(test_humanstress_exe clm csm_share)
