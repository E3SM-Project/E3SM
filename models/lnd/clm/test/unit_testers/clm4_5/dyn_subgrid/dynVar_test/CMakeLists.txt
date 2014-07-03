set (pfunit_sources
  test_dynVarTimeUninterp.pf
  test_dynVarTimeInterp.pf)

# extra sources used for this test, which are not .pf files
set (extra_sources
  test_dynVarShared.F90)

create_pfUnit_test(dynVar test_dynVar_exe
  "${pfunit_sources}" "${extra_sources}")

target_link_libraries(test_dynVar_exe clm csm_share)