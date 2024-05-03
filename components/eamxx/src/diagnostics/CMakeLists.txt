set(DIAGNOSTIC_SRCS
  atm_density.cpp
  dry_static_energy.cpp
  exner.cpp
  field_at_height.cpp
  field_at_level.cpp
  field_at_pressure_level.cpp
  longwave_cloud_forcing.cpp
  potential_temperature.cpp
  precip_surf_mass_flux.cpp
  relative_humidity.cpp
  sea_level_pressure.cpp
  shortwave_cloud_forcing.cpp
  surf_upward_latent_heat_flux.cpp
  vapor_flux.cpp
  vertical_layer.cpp
  virtual_temperature.cpp
  water_path.cpp
  wind_speed.cpp
  aodvis.cpp
  number_path.cpp
  aerocom_cld.cpp
  atm_tend.cpp
)

add_library(diagnostics ${DIAGNOSTIC_SRCS})
target_link_libraries(diagnostics PUBLIC scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
