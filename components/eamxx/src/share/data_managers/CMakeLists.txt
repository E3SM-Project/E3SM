add_library(eamxx_data_managers
  field_manager.cpp
  grids_manager.cpp
  IOPDataManager.cpp
)

target_link_libraries(eamxx_data_managers PUBLIC
  eamxx_field
  eamxx_grid
  eamxx_remap
  eamxx_scorpio_interface
)

if (NOT SCREAM_LIB_ONLY)
  # Add mesh free grids manager to the data managers lib (requires linking prop checks too)
  target_sources(eamxx_data_managers PRIVATE mesh_free_grids_manager.cpp)
  target_link_libraries(eamxx_data_managers PUBLIC
    eamxx_property_checks)

  add_subdirectory(tests)
endif()
