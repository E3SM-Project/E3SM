set(DIAGNOSTIC_SRCS
  aerocom_cld.cpp
  aodvis.cpp
  atm_backtend.cpp
  atm_density.cpp
  dry_static_energy.cpp
  exner.cpp
  field_at_height.cpp
  field_at_level.cpp
  field_at_pressure_level.cpp
  horiz_avg.cpp
  longwave_cloud_forcing.cpp
  number_path.cpp
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
)

add_library(diagnostics ${DIAGNOSTIC_SRCS})
target_link_libraries(diagnostics PUBLIC scream_share)

if (NOT SCREAM_LIB_ONLY AND NOT SCREAM_ONLY_GENERATE_BASELINES)
  add_subdirectory(tests)
endif()
