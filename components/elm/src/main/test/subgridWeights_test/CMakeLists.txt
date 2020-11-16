create_pFUnit_test(subgridWeights test_subgridWeights_exe
  "test_subgridWeights.pf" "")

target_link_libraries(test_subgridWeights_exe clm csm_share)