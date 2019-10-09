# Local pFUnit files.
set(pf_sources
  test_shr_cal.pf)

# Sources to test.
set(sources_needed
  shr_cal_mod.F90
  shr_kind_mod.F90
  shr_const_mod.F90
  shr_sys_mod.nompi_abortthrows.F90
  shr_abort_mod.abortthrows.F90
  shr_string_mod.F90
  shr_strconvert_mod.F90
  shr_infnan_mod.F90
  shr_log_mod.F90
  shr_assert_mod.F90
  shr_timer_mod.F90)

extract_sources("${sources_needed}" "${share_sources}" test_sources)

# shr_cal_mod depends on the esmf wrf time manager. So we need to add
# those sources to the list, too. (We can't build those into a separate
# library because they depend on shr_sys_mod - i.e., there are circular
# dependencies.)
list(APPEND test_sources ${esmf_wrf_timemgr_sources})

# This is needed to make ESMF_Stubs.F90 avoid pulling in mpi
add_definitions(-DHIDE_MPI)

# Do source preprocessing and add the executable.
create_pFUnit_test(shr_cal_mod shr_cal_exe "${pf_sources}"
  "${test_sources}")

declare_generated_dependencies(shr_cal_exe "${share_genf90_sources}")