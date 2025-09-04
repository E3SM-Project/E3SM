add_library(eamxx_grid
  grid/abstract_grid.cpp
  grid/grids_manager.cpp
  grid/grid_import_export.cpp
  grid/se_grid.cpp
  grid/point_grid.cpp
  grid/mesh_free_grids_manager.cpp
)

target_link_libraries(eamxx_grid PUBLIC eamxx_field)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
