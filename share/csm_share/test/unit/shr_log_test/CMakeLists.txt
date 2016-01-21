set(pf_sources
  test_error_printers.pf)

set(sources_needed shr_kind_mod.F90 shr_infnan_mod.F90
  shr_strconvert_mod.F90 shr_log_mod.F90)

extract_sources("${sources_needed}" "${share_sources}" test_sources)

create_pFUnit_test(shr_log shr_log_exe "${pf_sources}"
  "${test_sources}")

declare_generated_dependencies(shr_log_exe "${share_genf90_sources}")
