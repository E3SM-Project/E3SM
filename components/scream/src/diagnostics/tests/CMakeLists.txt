# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  include(ScreamUtils)

  set( NEED_LIBS scream_share diagnostics )

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

endif()
