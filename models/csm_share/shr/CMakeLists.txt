set(genf90_files shr_infnan_mod.F90.in shr_assert_mod.F90.in)

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_SOURCE_DIR}
  share_genf90_sources)

sourcelist_to_parent(share_genf90_sources)

list(APPEND share_sources "${share_genf90_sources}")

list(APPEND share_sources shr_kind_mod.F90 shr_const_mod.F90
  shr_sys_mod.F90 shr_log_mod.F90 shr_spfn_mod.F90)

sourcelist_to_parent(share_sources)
