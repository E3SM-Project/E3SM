include (ScreamUtils)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D
CreateDynamicsLib("theta-l_kokkos"  4   72   35)

# Create the test
SET (TEST_LABELS "dynamics;driver;shoc;cld;p3;rrtmgp;physics")
set (NEED_LIBS cld_fraction shoc spa p3 scream_rrtmgp ${dynLibName} scream_control scream_share physics_share)
CreateUnitTest(homme_shoc_cld_spa_p3_rrtmgp "homme_shoc_cld_spa_p3_rrtmgp.cpp" "${NEED_LIBS}"
               LABELS ${TEST_LABELS} 
               MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
               PROPERTIES FIXTURES_SETUP homme_shoc_cld_spa_p3_rrtmgp_generate_output_nc_files
)
# Make sure the SPA data directory is correctly defined
if (NOT EXISTS ${SCREAM_SPA_DATA_DIR})
  message(FATAL_ERROR "Failed, SCREAM_SPA_DATA_DIR = ${SCREAM_SPA_DATA_DIR} does not exist")
endif()
if (NOT IS_DIRECTORY ${SCREAM_SPA_DATA_DIR})
  message(FATAL_ERROR "Failed, SCREAM_SPA_DATA_DIR = ${SCREAM_SPA_DATA_DIR} is not a directory")
endif() 
if (NOT IS_ABSOLUTE ${SCREAM_SPA_DATA_DIR})
  message(FATAL_ERROR "Failed, SCREAM_SPA_DATA_DIR = ${SCREAM_SPA_DATA_DIR} is not an absolute directory")
endif() 

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestProfile(NUM_STEPS 2 4 48)  # 1h 2h 24h

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/homme_shoc_cld_spa_p3_rrtmgp_output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/homme_shoc_cld_spa_p3_rrtmgp_output.yaml)

# Set homme's test options, so that we can configure the namelist correctly
# Discretization/algorithm settings
set (HOMME_TEST_NE 2)
set (HOMME_TEST_LIM 9)
set (HOMME_TEST_REMAP_FACTOR 3)
set (HOMME_TEST_TRACERS_FACTOR 1)
set (HOMME_TEST_TIME_STEP 300)
SET (HOMME_THETA_FORM 1)
SET (HOMME_TTYPE 5)
SET (HOMME_SE_FTYPE 0)

# Hyperviscosity settings
set (HOMME_TEST_HVSCALING 0)
SET (HOMME_TEST_HVS 1)
SET (HOMME_TEST_HVS_TOM 0)

set (HOMME_TEST_NU 7e15)
set (HOMME_TEST_NUDIV 1e15)
set (HOMME_TEST_NUTOP 2.5e5)

# Testcase settings
set (HOMME_TEST_MOISTURE notdry)
SET (HOMME_THETA_HY_MODE true)

# Vert coord settings
set (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
set (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)

# Configure the namelist into the test directory
configure_file(${SCREAM_SOURCE_DIR}/src/dynamics/homme/tests/theta.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Copy initial conditions files
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/homme_shoc_cld_spa_p3_rrtmgp_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/homme_shoc_cld_spa_p3_rrtmgp_init_ne2np4.nc COPYONLY)

# Create vcoord directory
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vcoord)

# Copy vcoord files to build dir
configure_file(${HOMME_SOURCE_DIR}/test/vcoord/acme-72i.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
configure_file(${HOMME_SOURCE_DIR}/test/vcoord/acme-72m.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)

# P3 data files
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data/p3_lookup_table_1.dat-v4.1.1 COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data/p3_lookup_table_2.dat-v4.1.1 COPYONLY)
# RRTMGP data files
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)

## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if max mpi ranks is >1
# This test requires CPRNC
if (SCREAM_TEST_MAX_RANKS GREATER 1)
  include (BuildCprnc)
  BuildCprnc()
  SET (BASE_TEST_NAME "homme_shoc_cld_spa_p3_rrtmgp")
  foreach (MPI_RANKS RANGE 2 ${SCREAM_TEST_MAX_RANKS})
    set (SRC_FILE "${BASE_TEST_NAME}_output.INSTANT.Steps_x${NUM_STEPS}.np${MPI_RANKS}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_output.INSTANT.Steps_x${NUM_STEPS}.np1.nc")
    configure_file (${SCREAM_BASE_DIR}/cmake/CprncTest.cmake
                    ${CMAKE_CURRENT_BINARY_DIR}/CprncTest_np${MPI_RANKS}.cmake @ONLY)
    set(TEST_NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME} 
              COMMAND cmake -P CprncTest_np${MPI_RANKS}.cmake
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED homme_shoc_cld_spa_p3_rrtmgp_generate_output_nc_files)
  endforeach()
endif()
