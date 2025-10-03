add_library(eamxx_remap
  abstract_remapper.cpp
  coarsening_remapper.cpp
  horiz_interp_remapper_base.cpp
  horiz_interp_remapper_data.cpp
  iop_remapper.cpp
  identity_remapper.cpp
  refining_remapper_p2p.cpp
  vertical_remapper.cpp
)

if (EAMXX_ENABLE_EXPERIMENTAL_CODE)
  target_sources(eamxx_remap PRIVATE
    refining_remapper_rma.cpp
  )
endif()

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
