set(sources_needed circle.F90)
extract_sources("${sources_needed}" "${circle_area_sources}" test_sources)

create_pFUnit_test(pFunit_circle_area pFunittest_circle_area_exe "test_circle.pf" ${test_sources})

