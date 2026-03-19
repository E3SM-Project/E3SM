add_library(eamxx_remap
  abstract_remapper.cpp
  horiz_interp_remapper_data.cpp
  horizontal_remapper.cpp
  iop_remapper.cpp
  identity_remapper.cpp
  vertical_remapper.cpp
)

target_link_libraries(eamxx_remap PUBLIC
  eamxx_core
  eamxx_scorpio_interface
  eamxx_utils
  eamxx_field
  eamxx_grid
  ekat::Algorithm)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
