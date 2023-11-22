include (ScreamUtils)

set (BASE_TEST_NAME homme_shoc_cld_spa_p3_rrtmgp_128levels)
set (FIXTURES_BASE_NAME ${BASE_TEST_NAME}_generate_output_nc_files)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D
CreateDynamicsLib("theta-l_kokkos"  4   128   10)

# Create the test
CreateADUnitTest(${BASE_TEST_NAME}
  LABELS dynamics shoc cld p3 rrtmgp physics
  LIBS cld_fraction shoc nudging spa p3 scream_rrtmgp ${dynLibName} diagnostics
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
set (HOMME_TEST_VCOORD_INT_FILE acme-128i.ascii)
set (HOMME_TEST_VCOORD_MID_FILE acme-128m.ascii)

# Configure the namelist into the test directory
configure_file(${SCREAM_SRC_DIR}/dynamics/homme/tests/theta.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  scream/init/spa_file_unified_and_complete_ne4_20220428.nc
  scream/init/${EAMxx_tests_IC_FILE_128lev}
  cam/topo/${EAMxx_tests_TOPO_FILE}
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()


## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if running with 2+ ranks configurations
# This test requires CPRNC
set (BASE_TEST_NAME "homme_shoc_cld_spa_p3_rrtmgp_128levels")
set (FIXTURES_BASE_NAME homme_shoc_cld_spa_p3_rrtmgp_128levels_generate_output_nc_files)
if (TEST_RANK_END GREATER TEST_RANK_START)
  include (BuildCprnc)
  BuildCprnc()

  math (EXPR CMP_RANK_START ${TEST_RANK_START}+1)
  foreach (MPI_RANKS RANGE ${CMP_RANK_START} ${TEST_RANK_END})

    set (SRC_FILE "${BASE_TEST_NAME}_output.INSTANT.nsteps_x1.np${MPI_RANKS}.${RUN_T0}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_output.INSTANT.nsteps_x1.np${TEST_RANK_START}.${RUN_T0}.nc")
    set (TEST_NAME "${BASE_TEST_NAME}_np${TEST_RANK_START}_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES
        LABELS "dynamics;shoc;cld;p3;rrtmgp;physics"
        FIXTURES_REQUIRED "${FIXTURES_BASE_NAME}_np${MPI_RANKS}_omp1;${FIXTURES_BASE_NAME}_np1_omp1")
  endforeach()
endif()

# Check tendency calculation
foreach (NRANKS RANGE ${TEST_RANK_START} ${TEST_RANK_END})
  set (script ${SCREAM_BASE_DIR}/scripts/check-tendencies)
  set (fname ${BASE_TEST_NAME}_output.INSTANT.nsteps_x1.np${NRANKS}.${RUN_T0}.nc)
  set (tname ${BASE_TEST_NAME}_tend_check_np${NRANKS})
  add_test (NAME ${tname}
    COMMAND ${script} -f ${fname} -v qv T_mid -t EAMxx_qv_tend EAMxx_T_mid_tend
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (${tname} PROPERTIES
      LABELS "dynamics;shoc;cld;p3;rrtmgp;physics"
      FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${NRANKS}_omp1)
endforeach()
