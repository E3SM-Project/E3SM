include (ScreamUtils)

set (TEST_BASE_NAME homme_shoc_cld_spa_p3_rrtmgp_dp)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D
CreateDynamicsLib("theta-l_kokkos"  4   128   10)

# Create the test
set (TEST_LABELS "dynamics;driver;shoc;cld;spa;p3;rrtmgp;physics;dp")
CreateUnitTest(homme_shoc_cld_spa_p3_rrtmgp_dp "homme_shoc_cld_spa_p3_rrtmgp_dp.cpp"
  LABELS ${TEST_LABELS}
  LIBS cld_fraction shoc spa p3 scream_rrtmgp ${dynLibName} scream_control diagnostics
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
set (ATM_TIME_STEP 100)
SetVarDependingOnTestSize(NUM_STEPS 2 4 48)  # 1h 2h 24h
set (RUN_T0 1999-07-10-00000)

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
set (HOMME_TEST_VCOORD_INT_FILE acme-128i.ascii)
set (HOMME_TEST_VCOORD_MID_FILE acme-128m.ascii)

# Configure the namelist into the test directory
configure_file(${SCREAM_SRC_DIR}/dynamics/homme/tests/theta-dp.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  scream/init/${EAMxx_tests_IC_FILE_128lev}
  cam/topo/${EAMxx_tests_TOPO_FILE}
  cam/scam/iop/DYCOMSrf01_iopfile_4scam.nc
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
  LABELS dynamics physics shoc cld spa p3 rrtmgp dp
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_END}.${RUN_T0}.nc)
  set (SRC_FILE ${CMAKE_CURRENT_BINARY_DIR}/${OUT_FILE})
  set (TGT_FILE ${SCREAM_BASELINES_DIR}/data/${OUT_FILE})
  add_test (
    NAME ${TEST_BASE_NAME}_baseline_cmp
    COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(
    ${TEST_BASE_NAME}_baseline_cmp
    PROPERTIES
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${TEST_RANK_END}_omp1)

  # The np${TEST_RANK_END} test is used for baselines tests, so add it to the baseline label
  if (TEST_RANK_END GREATER 1)
    # ECUT only adds _npX if we have more than one value of X, or if X>1
    set (TEST_SUFFIX _np${TEST_RANK_END})
  endif()
  set_tests_properties (${TEST_BASE_NAME}${TEST_SUFFIX} PROPERTIES LABELS baseline)

  # test-all-scream will read this file to get the list of baseline nc files to
  # copy into the baseline dir
  file (APPEND ${SCREAM_TEST_OUTPUT_DIR}/baseline_list
    "${SRC_FILE}\n"
  )
endif()
