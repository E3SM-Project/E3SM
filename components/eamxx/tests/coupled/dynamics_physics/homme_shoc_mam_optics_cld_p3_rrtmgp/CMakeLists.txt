include (ScreamUtils)

set (TEST_BASE_NAME homme_shoc_mam_optics_cld_p3_rrtmgp)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D
CreateDynamicsLib("theta-l_kokkos"  4   72   70)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LIBS cld_fraction ${dynLibName} shoc p3 scream_rrtmgp mam
  LABELS dynamics shoc cld p3 rrtmgp physics mam4_optics
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 4 48)  # 1h 2h 24h
set (RUN_T0 2021-10-12-45000)

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
set (HOMME_TEST_NE 2)
set (HOMME_TEST_LIM 9)
set (HOMME_TEST_REMAP_FACTOR 3)
set (HOMME_TEST_TRACERS_FACTOR 1)
set (HOMME_TEST_TIME_STEP 300)
set (HOMME_THETA_FORM 1)
set (HOMME_TTYPE 5)
set (HOMME_SE_FTYPE 0)
set (HOMME_TEST_TRANSPORT_ALG 0)
set (HOMME_TEST_CUBED_SPHERE_MAP 0)

# Hyperviscosity settings
set (HOMME_TEST_HVSCALING 0)
set (HOMME_TEST_HVS 1)
set (HOMME_TEST_HVS_TOM 0)
set (HOMME_TEST_HVS_Q 1)

set (HOMME_TEST_NU 7e15)
set (HOMME_TEST_NUDIV 1e15)
set (HOMME_TEST_NUTOP 2.5e5)

# Testcase settings
set (HOMME_TEST_MOISTURE notdry)
set (HOMME_THETA_HY_MODE true)

# Vert coord settings
set (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
set (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)

# Configure the namelist into the test directory
configure_file(${SCREAM_SRC_DIR}/dynamics/homme/tests/theta.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  # scream/init/spa_init_ne2np4.nc
  # scream/init/spa_file_unified_and_complete_ne2np4L72_20231222.nc
  scream/init/${EAMxx_tests_IC_FILE_72lev}
  cam/topo/${EAMxx_tests_TOPO_FILE}
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()
#GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
#GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x${NUM_STEPS}.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS dynamics physics shoc cld p3 rrtmgp mam4_optics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)
