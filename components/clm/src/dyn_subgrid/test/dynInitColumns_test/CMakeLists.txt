create_pFUnit_test(dynInitColumns test_dynInitColumns_exe
  "test_init_columns.pf" "")

target_link_libraries(test_dynInitColumns_exe clm csm_share)
