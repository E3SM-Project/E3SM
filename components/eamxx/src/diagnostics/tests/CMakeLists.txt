include(ScreamUtils)

function (createDiagTest test_name test_srcs)
  CreateUnitTest(${test_name} "${test_srcs}"
    LIBS diagnostics physics_share
    LABELS diagnostics)
endfunction ()

# Test extracting a single level of a field
CreateDiagTest(field_at_level "field_at_level_tests.cpp")

# Test interpolating a field onto a single pressure level
CreateDiagTest(field_at_pressure_level "field_at_pressure_level_tests.cpp")

# Test interpolating a field at a specific height
CreateDiagTest(field_at_height "field_at_height_tests.cpp")

# Test potential temperature diagnostic
CreateDiagTest(potential_temperature "potential_temperature_test.cpp")

# Test exner diagnostic
CreateDiagTest(exner_function "exner_test.cpp")

# Test virtual temperature
CreateDiagTest(virtual_temperature "virtual_temperature_test.cpp")

# Test atmosphere density
CreateDiagTest(atmosphere_density "atm_density_test.cpp")

# Test vertical layer (dz, z_int, z_mid)
CreateDiagTest(vertical_layer "vertical_layer_tests.cpp")

# Test dry static energy
CreateDiagTest(dry_static_energy "dry_static_energy_test.cpp")

# Test sea level pressure
CreateDiagTest(sea_level_pressure "sea_level_pressure_test.cpp")

# Test total water path
CreateDiagTest(water_path "water_path_tests.cpp")

# Test shortwave cloud forcing
CreateDiagTest(shortwave_cloud_forcing "shortwave_cloud_forcing_tests.cpp")

# Test longwave cloud forcing
CreateDiagTest(longwave_cloud_forcing "longwave_cloud_forcing_tests.cpp")

# Test Relative Humidity
CreateDiagTest(relative_humidity "relative_humidity_tests.cpp")

# Test Vapor Flux
CreateDiagTest(vapor_flux "vapor_flux_tests.cpp")

# Test precipitation mass surface flux
CreateDiagTest(precip_surf_mass_flux "precip_surf_mass_flux_tests.cpp")

# Test surface latent heat flux
CreateDiagTest(surface_upward_latent_heat_flux "surf_upward_latent_heat_flux_tests.cpp")

# Test wind speed diagnostic
CreateDiagTest(wind_speed "wind_speed_tests.cpp")

# Test AODVIS
CreateDiagTest(aodvis "aodvis_test.cpp")

# Test "number" paths
CreateDiagTest(number_paths "number_paths_tests.cpp")

# Test AEROCOM_CLD
CreateDiagTest(aerocom_cld "aerocom_cld_test.cpp")

# Test atm_tend
CreateDiagTest(atm_backtend "atm_backtend_test.cpp")

# Test horizontal averaging
CreateDiagTest(horiz_avg "horiz_avg_test.cpp")
