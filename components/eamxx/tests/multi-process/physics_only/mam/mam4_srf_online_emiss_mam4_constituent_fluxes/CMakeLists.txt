INCLUDE (ScreamUtils)

set (TEST_BASE_NAME mam4_srf_online_emiss_mam4_constituent_fluxes)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LIBS mam
  LABELS physics mam4_srf_online_emiss mam4_constituent_fluxes
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)  # 1h 4h 24h
set (RUN_T0 2021-10-12-45000)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_MAM4xx_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  scream/mam4xx/emissions/ne2np4/surface/DMSflux.2010.ne2np4_conserv.POPmonthlyClimFromACES4BGC_c20240726.nc
  scream/mam4xx/emissions/ne2np4/surface/cmip6_mam4_so2_surf_ne2np4_2010_clim_c20240723.nc
  scream/mam4xx/emissions/ne2np4/surface/cmip6_mam4_bc_a4_surf_ne2np4_2010_clim_c20240726.nc
  scream/mam4xx/emissions/ne2np4/surface/cmip6_mam4_num_a1_surf_ne2np4_2010_clim_c20240726.nc
  scream/mam4xx/emissions/ne2np4/surface/cmip6_mam4_num_a2_surf_ne2np4_2010_clim_c20240726.nc
  scream/mam4xx/emissions/ne2np4/surface/cmip6_mam4_num_a4_surf_ne2np4_2010_clim_c20240726.nc
  scream/mam4xx/emissions/ne2np4/surface/cmip6_mam4_pom_a4_surf_ne2np4_2010_clim_c20240726.nc
  scream/mam4xx/emissions/ne2np4/surface/cmip6_mam4_so4_a1_surf_ne2np4_2010_clim_c20240726.nc
  scream/mam4xx/emissions/ne2np4/surface/cmip6_mam4_so4_a2_surf_ne2np4_2010_clim_c20240726.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)
CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x${NUM_STEPS}.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS physics mam4_srf_online_emiss mam4_constituent_fluxes
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)
