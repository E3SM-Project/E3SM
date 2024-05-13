include (ScreamUtils)

set (TEST_BASE_NAME mam4_optics_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LABELS mam4_optics physics
  LIBS mam
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 12 24 36)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)


# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_MAM4xx_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

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

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x2.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS mam4_optics physics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x2.np${TEST_RANK_END}.${RUN_T0}.nc)
  CreateBaselineTest(${TEST_BASE_NAME} ${TEST_RANK_END} ${OUT_FILE} ${FIXTURES_BASE_NAME})
endif()