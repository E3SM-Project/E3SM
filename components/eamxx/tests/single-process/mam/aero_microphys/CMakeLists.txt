include (ScreamUtils)

set (TEST_BASE_NAME mam4_aero_microphys_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LABELS mam4_aero_microphys physics
  LIBS mam
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)  # 1h 2.5h 24h
set (RUN_T0 2021-10-12-45000)

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
     scream/init/${EAMxx_tests_IC_FILE_MAM4xx_72lev}
     scream/mam4xx/linoz/ne2np4/linoz1850-2015_2010JPL_CMIP6_10deg_58km_ne2np4_c20240724.nc
     scream/mam4xx/invariants/ne2np4/oxid_ne2np4_L26_1850-2015_c20240827.nc
     scream/mam4xx/linoz/Linoz_Chlorine_Loading_CMIP6_0003-2017_c20171114.nc
     scream/mam4xx/photolysis/RSF_GT200nm_v3.0_c080811.nc
     scream/mam4xx/photolysis/temp_prs_GT200nm_JPL10_c130206.nc
     scream/mam4xx/emissions/ne2np4/elevated/cmip6_mam4_so2_elev_ne2np4_2010_clim_c20240726.nc
     scream/mam4xx/emissions/ne2np4/elevated/cmip6_mam4_so4_a1_elev_ne2np4_2010_clim_c20240823.nc
     scream/mam4xx/emissions/ne2np4/elevated/cmip6_mam4_so4_a2_elev_ne2np4_2010_clim_c20240823.nc
     scream/mam4xx/emissions/ne2np4/elevated/cmip6_mam4_pom_a4_elev_ne2np4_2010_clim_c20240823.nc
     scream/mam4xx/emissions/ne2np4/elevated/cmip6_mam4_bc_a4_elev_ne2np4_2010_clim_c20240823.nc
     scream/mam4xx/emissions/ne2np4/elevated/cmip6_mam4_num_a1_elev_ne2np4_2010_clim_c20240823.nc
     scream/mam4xx/emissions/ne2np4/elevated/cmip6_mam4_num_a2_elev_ne2np4_2010_clim_c20240823.nc
     scream/mam4xx/emissions/ne2np4/elevated/cmip6_mam4_num_a4_elev_ne2np4_2010_clim_c20240823.nc
     scream/mam4xx/emissions/ne2np4/elevated/cmip6_mam4_soag_elev_ne2np4_2010_clim_c20240823.nc
     scream/mam4xx/drydep/season_wes.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()


# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS mam4_aero_microphys physics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.np${TEST_RANK_END}.${RUN_T0}.nc)
  CreateBaselineTest(${TEST_BASE_NAME} ${TEST_RANK_END} ${OUT_FILE} ${FIXTURES_BASE_NAME})
endif()
