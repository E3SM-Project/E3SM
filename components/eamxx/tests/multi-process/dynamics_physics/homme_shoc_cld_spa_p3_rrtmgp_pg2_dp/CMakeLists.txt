include (ScreamUtils)

set (TEST_BASE_NAME homme_shoc_cld_spa_p3_rrtmgp_pg2_dp)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D
CreateDynamicsLib("theta-l_kokkos"  4   72   10)

# Create the test
set (TEST_LABELS "dynamics;driver;tms;shoc;cld;spa;p3;rrtmgp;physics;dp")
CreateUnitTest(homme_shoc_cld_spa_p3_rrtmgp_pg2_dp "homme_shoc_cld_spa_p3_rrtmgp_pg2_dp.cpp"
  LABELS ${TEST_LABELS}
  LIBS cld_fraction tms shoc spa iop_forcing p3 scream_rrtmgp ${dynLibName} scream_control diagnostics
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
set (ATM_TIME_STEP 100)
SetVarDependingOnTestSize(NUM_STEPS 2 4 48)  # 1h 2h 24h
set (RUN_T0 2000-07-24-61100)

# Determine num subcycles needed to keep shoc dt<=300s
set (SHOC_MAX_DT 300)
math (EXPR MAC_MIC_SUBCYCLES "(${ATM_TIME_STEP} + ${SHOC_MAX_DT} - 1) / ${SHOC_MAX_DT}")

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Set homme's test options, so that we can configure the namelist correctly
# Discretization/algorithm settings
set (HOMME_TEST_LIM 9)
set (HOMME_TEST_REMAP_FACTOR 2)
set (HOMME_TEST_TRACERS_FACTOR 1)
set (HOMME_TEST_TIME_STEP 8.3333333333333)
set (HOMME_THETA_FORM 1)
set (HOMME_TTYPE 9)
set (HOMME_SE_FTYPE 2)
set (HOMME_TEST_TRANSPORT_ALG 0)
set (HOMME_TEST_CUBED_SPHERE_MAP 2)
set (HOMME_TEST_NE_X 5)
set (HOMME_TEST_NE_Y 5)
set (HOMME_TEST_LX 50000)
set (HOMME_TEST_LY 50000)

# Hyperviscosity settings
set (HOMME_TEST_HVSCALING 3)
set (HOMME_TEST_HVS 1)
set (HOMME_TEST_HVS_TOM 1)
set (HOMME_TEST_HVS_Q 1)

set (HOMME_TEST_NU 0.216784)
set (HOMME_TEST_NUTOP 1e4)

# Testcase settings
set (HOMME_TEST_MOISTURE notdry)
set (HOMME_THETA_HY_MODE false)

# Vert coord settings
set (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
set (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)

# Configure the namelist into the test directory
configure_file(${SCREAM_SRC_DIR}/dynamics/homme/tests/theta-dp.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  scream/init/${EAMxx_tests_IC_FILE_72lev}
  cam/topo/${EAMxx_tests_TOPO_FILE}
  cam/scam/iop/CASS_iopfile_4scam.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x${NUM_STEPS}.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS dynamics physics tms shoc cld spa p3 rrtmgp dp pg2
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.${RUN_T0}.nc)
  CreateBaselineTest(${TEST_BASE_NAME} ${TEST_RANK_END} ${OUT_FILE} ${FIXTURES_BASE_NAME})
endif()
