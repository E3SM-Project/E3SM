set(genf90_files pio.F90.in)

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_BINARY_DIR}
  pio_genf90_sources)

sourcelist_to_parent(pio_genf90_sources)

list(APPEND pio_sources "${pio_genf90_sources}")

sourcelist_to_parent(pio_sources)