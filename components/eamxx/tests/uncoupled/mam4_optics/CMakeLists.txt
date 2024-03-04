include (ScreamUtils)

# This test needs to be reworked for microphysics -- currently it's still using
# input for nucleation.

# Create the test
CreateADUnitTest(mam4_optics_standalone
  LABELS mam4_optics physics
  LIBS mam
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
)

# Set AD configurable options
set (ATM_TIME_STEP 1)
SetVarDependingOnTestSize(NUM_STEPS 12 24 36)
set (RUN_T0 2021-10-12-45000)

# ## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  scream/mam4xx/physprops/mam4_mode1_rrtmg_aeronetdust_c20240206.nc
  scream/mam4xx/physprops/mam4_mode2_rrtmg_c20240206.nc
  scream/mam4xx/physprops/mam4_mode3_rrtmg_aeronetdust_c20240206.nc
  scream/mam4xx/physprops/mam4_mode4_rrtmg_c20240206.nc
  scream/mam4xx/physprops/water_refindex_rrtmg_c20240206.nc
  scream/mam4xx/physprops/ocphi_rrtmg_c20240206.nc
  scream/mam4xx/physprops/dust_aeronet_rrtmg_c20240206.nc
  scream/mam4xx/physprops/ssam_rrtmg_c20240206.nc
  scream/mam4xx/physprops/sulfate_rrtmg_c20240206.nc
  scream/mam4xx/physprops/ocpho_rrtmg_c20240206.nc
  scream/mam4xx/physprops/bcpho_rrtmg_c20240206.nc
  scream/mam4xx/physprops/poly_rrtmg_c20240206.nc
  scream/init/scream_unit_tests_aerosol_optics_ne2np4L72_20220822.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()
