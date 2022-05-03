include (ScreamUtils)

# Create the test
set (TEST_LABELS "shoc;cld;spa;p3;rrtmgp;physics")
set (NEED_LIBS shoc cld_fraction spa p3 scream_rrtmgp scream_control scream_share physics_share)
CreateUnitTest(shoc_cld_spa_p3_rrtmgp shoc_cld_spa_p3_rrtmgp.cpp "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  PROPERTIES FIXTURES_SETUP shoc_cld_spa_p3_rrtmgp_generate_output_nc_files
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestProfile(NUM_STEPS 2 5 48)  # 1h 4h 24h

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  init/shoc_cld_spa_p3_rrtmgp_init_ne2np4.nc
  init/spa_init_ne2np4.nc
  init/map_ne4np4_to_ne2np4_mono.nc
  init/spa_file_unified_and_complete_ne4_20220428.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shoc_cld_spa_p3_rrtmgp_output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_cld_spa_p3_rrtmgp_output.yaml)

## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if max mpi ranks is >1
# This test requires CPRNC
if (SCREAM_TEST_MAX_RANKS GREATER 1)
  include (BuildCprnc)
  BuildCprnc()

  set (BASE_TEST_NAME "shoc_cld_spa_p3_rrtmgp")
  foreach (MPI_RANKS RANGE 2 ${SCREAM_TEST_MAX_RANKS})
    set (SRC_FILE "${BASE_TEST_NAME}_output.INSTANT.Steps_x${NUM_STEPS}.np${MPI_RANKS}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_output.INSTANT.Steps_x${NUM_STEPS}.np1.nc")
    set (TEST_NAME "${BASE_TEST_NAME}_np1_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED shoc_cld_spa_p3_rrtmgp_generate_output_nc_files)
  endforeach()
endif()
