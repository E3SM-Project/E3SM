add_library(eamxx_grid
  abstract_grid.cpp
  grid_import_export.cpp
  se_grid.cpp
  point_grid.cpp
)

target_link_libraries (eamxx_grid PUBLIC
  eamxx_core
  eamxx_utils
  eamxx_field)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
