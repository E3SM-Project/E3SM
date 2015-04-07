set(genf90_files
  ncdio_pio_fake.F90.in
  ncdio_var.F90.in
  )

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_SOURCE_DIR} clm_genf90_sources)

sourcelist_to_parent(clm_genf90_sources)

list(APPEND clm_sources "${clm_genf90_sources}")

list(APPEND clm_sources
  GetGlobalValuesMod_stub.F90
  histFileMod_stub.F90
  ncdio_pio_fake.F90
  ncdio_var.F90
  )

sourcelist_to_parent(clm_sources)
