create_pFUnit_test(avect_wrapper avect_wrapper_exe
  "test_avect_wrapper.pf" "")

target_link_libraries(avect_wrapper_exe ${DRV_UNIT_TEST_LIBS})
