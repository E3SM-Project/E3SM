set(pf_sources
  test_toString.pf)

set(sources_needed shr_kind_mod.F90 shr_infnan_mod.F90
  shr_strconvert_mod.F90)

extract_sources("${sources_needed}" "${share_sources}" test_sources)

create_pFUnit_test(shr_strconvert shr_strconvert_exe "${pf_sources}"
  "${test_sources}")

declare_generated_dependencies(shr_strconvert_exe "${share_genf90_sources}")
