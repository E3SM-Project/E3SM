# Local pFUnit files.
set(pf_sources
  test_tridiagonal_inverse.pf)

# Sources to test.
set(sources_needed
  shr_assert_mod.F90
  shr_infnan_mod.F90
  shr_kind_mod.F90
  shr_log_mod.F90
  shr_matrix_mod.F90
  shr_strconvert_mod.F90
  shr_sys_mod.nompi_abortthrows.F90)
extract_sources("${sources_needed}" "${share_sources}" test_sources)

# Do source preprocessing and add the executable.
create_pFUnit_test(shr_matrix_mod shr_matrix_mod_exe
  "${pf_sources}" "${test_sources}")

declare_generated_dependencies(shr_matrix_mod_exe "${share_genf90_sources}")