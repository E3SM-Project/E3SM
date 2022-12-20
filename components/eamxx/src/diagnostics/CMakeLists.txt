set(DIAGNOSTIC_SRCS
  atm_density.cpp
  dry_static_energy.cpp
  exner.cpp
  field_at_level.cpp
  field_at_pressure_level.cpp
  ice_cloud_mask.cpp
  ice_water_path.cpp
  liquid_water_path.cpp
  longwave_cloud_forcing.cpp
  meridional_vapor_flux.cpp
  potential_temperature.cpp
  precip_ice_surf_mass_flux.cpp
  precip_liq_surf_mass_flux.cpp
  rain_water_path.cpp
  relative_humidity.cpp
  rime_water_path.cpp
  sea_level_pressure.cpp
  shortwave_cloud_forcing.cpp
  vapor_water_path.cpp
  vertical_layer_interface.cpp
  vertical_layer_midpoint.cpp
  vertical_layer_thickness.cpp
  virtual_temperature.cpp
  zonal_vapor_flux.cpp
)

add_library(diagnostics ${DIAGNOSTIC_SRCS})
target_include_directories(diagnostics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(diagnostics PUBLIC scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
