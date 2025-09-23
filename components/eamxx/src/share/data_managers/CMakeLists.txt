add_library(eamxx_managers
  field_manager.cpp
  grids_manager.cpp
  IOPDataManager.cpp
)

target_link_libraries(eamxx_managers PUBLIC
  eamxx_field
  eamxx_grid
  eamxx_remap
  eamxx_scorpio_interface
)

if (NOT SCREAM_LIB_ONLY)
  # It's used by quite a few unit tests
  target_sources(eamxx_managers PRIVATE mesh_free_grids_manager.cpp)

  add_subdirectory(tests)
endif()
