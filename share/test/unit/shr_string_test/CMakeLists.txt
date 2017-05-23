set (pf_sources
  test_shr_string.pf
  )

set(sources_needed
  shr_string_mod.F90
  shr_assert_mod.F90
  shr_infnan_mod.F90
  shr_kind_mod.F90
  shr_log_mod.F90
  shr_strconvert_mod.F90
  shr_sys_mod.nompi_abortthrows.F90
  shr_timer_mod.F90)

extract_sources("${sources_needed}" "${share_sources}" test_sources)

create_pFUnit_test(shr_string shr_string_exe "${pf_sources}" "${test_sources}")

declare_generated_dependencies(shr_string_exe "${share_genf90_sources}")