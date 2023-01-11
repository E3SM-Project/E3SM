include (ScreamUtils)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D
CreateDynamicsLib("theta-l_kokkos"  4   128   10)

# Create the test
set (TEST_LABELS "dynamics;driver;shoc;cld;p3;rrtmgp;physics")
set (NEED_LIBS cld_fraction shoc spa p3 scream_rrtmgp ${dynLibName} scream_control scream_share physics_share diagnostics)
CreateUnitTest(homme_shoc_cld_spa_p3_rrtmgp_128levels "homme_shoc_cld_spa_p3_rrtmgp_128levels.cpp" "${NEED_LIBS}"
  LABELS ${TEST_LABELS}
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  PROPERTIES FIXTURES_SETUP homme_shoc_cld_spa_p3_rrtmgp_128levels_generate_output_nc_files
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 4 48)  # 1h 2h 24h

# Determine num subcycles needed to keep shoc dt<=300s
set (SHOC_MAX_DT 300)
math (EXPR MAC_MIC_SUBCYCLES "(${ATM_TIME_STEP} + ${SHOC_MAX_DT} - 1) / ${SHOC_MAX_DT}")

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
             ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/homme_shoc_cld_spa_p3_rrtmgp_128levels_output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/homme_shoc_cld_spa_p3_rrtmgp_128levels_output.yaml)

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

# Hyperviscosity settings
set (HOMME_TEST_HVSCALING 0)
set (HOMME_TEST_HVS 1)
set (HOMME_TEST_HVS_TOM 0)

set (HOMME_TEST_NU 7e15)
set (HOMME_TEST_NUDIV 1e15)
set (HOMME_TEST_NUTOP 2.5e5)

# Testcase settings
set (HOMME_TEST_MOISTURE notdry)
set (HOMME_THETA_HY_MODE true)

# Vert coord settings
set (HOMME_TEST_VCOORD_INT_FILE acme-128i.ascii)
set (HOMME_TEST_VCOORD_MID_FILE acme-128m.ascii)

# Configure the namelist into the test directory
configure_file(${SCREAM_SRC_DIR}/dynamics/homme/tests/theta.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  scream/init/spa_file_unified_and_complete_ne4_20220428.nc
  scream/init/${EAMxx_tests_IC_FILE_128lev}
  cam/topo/USGS-gtopo30_ne4np4_16x.c20160612.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()


## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if running with 2+ ranks configurations
# This test requires CPRNC
if (TEST_RANK_END GREATER TEST_RANK_START)
  include (BuildCprnc)
  BuildCprnc()

  set (BASE_TEST_NAME "homme_shoc_cld_spa_p3_rrtmgp_128levels")
  foreach (MPI_RANKS RANGE ${TEST_RANK_START} ${TEST_RANK_END})
    set (SRC_FILE "${BASE_TEST_NAME}_output.INSTANT.nsteps_x${NUM_STEPS}.np${MPI_RANKS}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_output.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_START}.nc")
    set (TEST_NAME "${BASE_TEST_NAME}_np${TEST_RANK_START}_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED homme_shoc_cld_spa_p3_rrtmgp_128levels_generate_output_nc_files)
  endforeach()
endif()
