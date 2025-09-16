add_library(eamxx_data_managers
  field_manager.cpp
  grids_manager.cpp
)

target_link_libraries(eamxx_data_managers PUBLIC
  eamxx_field
  eamxx_grid
  eamxx_remap
)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
