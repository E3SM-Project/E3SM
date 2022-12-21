# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  include(ScreamUtils)

  set( NEED_LIBS scream_share diagnostics physics_share )

  # Test extracting a single level of a field
  CreateUnitTest(field_at_level "field_at_level_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )

  # Test interpolating a field onto a single pressure level
  CreateUnitTest(field_at_pressure_level "field_at_pressure_level_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )

  # Test potential temperature diagnostic
  CreateUnitTest(potential_temperature "potential_temperature_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test exner diagnostic
  CreateUnitTest(exner_function "exner_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test virtual temperature
  CreateUnitTest(virtual_temperature "virtual_temperature_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test atmosphere density
  CreateUnitTest(atmosphere_density "atm_density_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test vertical layer thickness (dz) 
  CreateUnitTest(vertical_layer_thickness "vertical_layer_thickness_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test vertical layer interface (z_int) 
  CreateUnitTest(vertical_layer_interface "vertical_layer_interface_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test vertical layer interface (z_mid) 
  CreateUnitTest(vertical_layer_midpoint "vertical_layer_midpoint_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test dry static energy 
  CreateUnitTest(dry_static_energy "dry_static_energy_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test sea level pressure 
  CreateUnitTest(sea_level_pressure "sea_level_pressure_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test total water path
  CreateUnitTest(water_path "water_path_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test shortwave cloud forcing
  CreateUnitTest(shortwave_cloud_forcing "shortwave_cloud_forcing_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test longwave cloud forcing
  CreateUnitTest(longwave_cloud_forcing "longwave_cloud_forcing_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test Relative Humidity 
  CreateUnitTest(relative_humidity "relative_humidity_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test Zonal Vapor Flux
  CreateUnitTest(zonal_vapor_flux "zonal_vapor_flux_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test Meridional Vapor Flux
  CreateUnitTest(meridional_vapor_flux "meridional_vapor_flux_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test Meridional Vapor Flux
  CreateUnitTest(ice_cloud_mask "ice_cloud_mask_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test precipitation ice mass surface flux
  CreateUnitTest(precip_ice_surf_mass_flux "precip_ice_surf_mass_flux_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )
  # Test precipitation liq mass surface flux
  CreateUnitTest(precip_liq_surf_mass_flux "precip_liq_surf_mass_flux_tests.cpp" "${NEED_LIBS}" LABELS "diagnostics" )

endif()
