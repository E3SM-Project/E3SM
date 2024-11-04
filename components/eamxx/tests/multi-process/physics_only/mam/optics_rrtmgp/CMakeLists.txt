INCLUDE (ScreamUtils)

set (TEST_BASE_NAME mam_optics_rrtmgp)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LIBS shoc scream_rrtmgp mam
  LABELS rrtmgp physics diagnostics mam4_optics
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)  # 1h 4h 24h
set (RUN_T0 2021-10-12-45000)

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
)

foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_MAM4xx_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

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
  LABELS rrtmgp mam4_optics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)
