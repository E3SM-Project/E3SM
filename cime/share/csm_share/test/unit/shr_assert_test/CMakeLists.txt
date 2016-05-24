set(pf_sources test_assert.pf test_assert_array.pf test_macro.pf
  test_ndebug.pf)

set(sources_needed shr_kind_mod.F90 shr_infnan_mod.F90
  shr_strconvert_mod.F90 shr_log_mod.F90
  shr_sys_mod.nompi_abortthrows.F90 shr_assert_mod.F90)

extract_sources("${sources_needed}" "${share_sources}" test_sources)

create_pFUnit_test(assert assert_test_exe "${pf_sources}" "${test_sources}")

declare_generated_dependencies(assert_test_exe "${share_genf90_sources}")

