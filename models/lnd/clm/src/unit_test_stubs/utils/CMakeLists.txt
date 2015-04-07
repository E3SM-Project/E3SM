set(genf90_files
  restUtilMod_stub.F90.in
  )

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_SOURCE_DIR} clm_genf90_sources)

sourcelist_to_parent(clm_genf90_sources)

list(APPEND clm_sources "${clm_genf90_sources}")

list(APPEND clm_sources
  clm_time_manager_stub.F90
  restUtilMod_stub.F90
  spmdMod_stub.F90
  )

sourcelist_to_parent(clm_sources)
