add_library(eamxx_property_checks
  property_check.cpp
  field_nan_check.cpp
  field_within_interval_check.cpp
  mass_and_energy_conservation_check.cpp
)

target_link_libraries(eamxx_property_checks PUBLIC
  eamxx_field
  eamxx_grid
  eamxx_physics_share
  ekat::Algorithm)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
