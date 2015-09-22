# Local pFUnit files.
set(pf_sources
  test_wv_sat.pf test_wv_sat_each_method.pf)

# Sources to test.
set(sources_needed
  shr_kind_mod.F90 shr_const_mod.F90 shr_wv_sat_mod.F90)
extract_sources("${sources_needed}" "${share_sources}" test_sources)

# Do source preprocessing and add the executable.
create_pFUnit_test(shr_wv_sat_mod shr_wv_sat_mod_exe "${pf_sources}"
  "${test_sources}")
