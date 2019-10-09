create_pFUnit_test(dynTimeInfo test_dynTimeInfo_exe
  "test_dynTimeInfo.pf" "")

target_link_libraries(test_dynTimeInfo_exe clm csm_share)