create_pFUnit_test(clm_glclnd test_clm_glclnd_exe
  "test_clm_glclnd.pf" "")

target_link_libraries(test_clm_glclnd_exe clm csm_share)