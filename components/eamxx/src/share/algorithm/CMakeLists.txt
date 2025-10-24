add_library(eamxx_algorithm
  eamxx_data_interpolation.cpp
  eamxx_time_interpolation.cpp
  eamxx_fv_phys_rrtmgp_active_gases_workaround.cpp
)

target_link_libraries(eamxx_algorithm PUBLIC
  eamxx_core
  eamxx_utils
  eamxx_field
  eamxx_remap
  eamxx_data_managers
  eamxx_io)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
