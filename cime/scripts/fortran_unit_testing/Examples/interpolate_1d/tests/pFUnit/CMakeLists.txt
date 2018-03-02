extract_sources("interpolate_1d.F90" "${interpolate_sources}" test_sources)
create_pFUnit_test(pFunit_interpolate pFunit_interpolate_exe "test_interpolate_point.pf" ${test_sources})
