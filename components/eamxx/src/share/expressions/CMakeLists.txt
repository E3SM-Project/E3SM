add_library(eamxx_expressions INTERFACE)

target_link_libraries (eamxx_expressions INTERFACE
  eamxx_field
  eamxx_utils
  ekat::KokkosUtils)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
