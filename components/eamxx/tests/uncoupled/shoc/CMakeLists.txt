include (ScreamUtils)

# Create the test
SET (TEST_LABELS "shoc;physics;driver")
set (NEED_LIBS shoc scream_control scream_share diagnostics)
CreateUnitTest(shoc_standalone "shoc_standalone.cpp" "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  PROPERTIES FIXTURES_SETUP shoc_generate_output_nc_files
)

# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Determine num subcycles needed to keep shoc dt<=300s
set (SHOC_MAX_DT 300)
math (EXPR NUM_SUBCYCLES "(${ATM_TIME_STEP} + ${SHOC_MAX_DT} - 1) / ${SHOC_MAX_DT}")

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/USGS-gtopo30_ne4np4pg2_16x_converted.c20200527.nc)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(shoc_standalone_output.yaml shoc_standalone_output.yaml)

## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if running with 2+ ranks configurations
# This test requires CPRNC
if (TEST_RANK_END GREATER TEST_RANK_START)
  include (BuildCprnc)
  BuildCprnc()
  SET (BASE_TEST_NAME "shoc")
  math (EXPR CMP_RANK_START ${TEST_RANK_START}+1)
  foreach (MPI_RANKS RANGE ${CMP_RANK_START} ${TEST_RANK_END})

    set (SRC_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.nsteps_x1.np${MPI_RANKS}.${RUN_T0}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_standalone_output.INSTANT.nsteps_x1.np${TEST_RANK_START}.${RUN_T0}.nc")

    set (TEST_NAME "${BASE_TEST_NAME}_np${TEST_RANK_START}_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED shoc_generate_output_nc_files)
  endforeach()
endif()

# Check tendency calculation
foreach (NRANKS RANGE ${TEST_RANK_START} ${TEST_RANK_END})
  set (script ${SCREAM_BASE_DIR}/scripts/check-tendencies)
  set (fname shoc_standalone_output.INSTANT.nsteps_x1.np${NRANKS}.${RUN_T0}.nc)
  set (tname shoc_tend_check_np${NRANKS})
  add_test (NAME ${tname}
    COMMAND ${script} -f ${fname}
    -v T_mid qv tke horiz_winds sgs_buoy_flux eddy_diff_mom qc cldfrac_liq
    -t shoc_T_mid_tend shoc_qv_tend shoc_tke_tend shoc_horiz_winds_tend shoc_sgs_buoy_flux_tend shoc_eddy_diff_mom_tend shoc_qc_tend shoc_cldfrac_liq_tend
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (${tname} PROPERTIES
    FIXTURES_REQUIRED shoc_generate_output_nc_files
    LABELS "${TEST_LABELS}")
endforeach()

# Check that the content of sliced vars U/V and corresponding surf_mom_flux
# components in all output NC files matches the content of corresponding components
# of the vector quantities
foreach (RANK RANGE ${TEST_RANK_START} ${TEST_RANK_END})
  set (nc_file shoc_standalone_output.INSTANT.nsteps_x1.np${RANK}.${RUN_T0}.nc)
  add_test (NAME check_U_V_slices_np${RANK}
            COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
            -s ${nc_file} -c "horiz_winds(:,:,1,:)=U(:,:,:)" "horiz_winds(:,:,2,:)=V(:,:,:)"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (check_U_V_slices_np${RANK} PROPERTIES
            FIXTURES_REQUIRED shoc_generate_output_nc_files)
  add_test (NAME check_surf_mom_flux_slices_np${RANK}
            COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
            -s ${nc_file} -c "surf_mom_flux(:,:,1)=surf_mom_flux_U(:,:)"
                             "surf_mom_flux(:,:,2)=surf_mom_flux_V(:,:)"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (check_surf_mom_flux_slices_np${RANK} PROPERTIES
            FIXTURES_REQUIRED shoc_generate_output_nc_files)
endforeach()

################################
#        MUST FAIL tests       #
################################

# These tests ensure that the compare-nc-files script doesn't always pass
set (nc_file shoc_standalone_output.INSTANT.nsteps_x1.np${TEST_RANK_START}.${RUN_T0}.nc)

# Legitimate diff
add_test (NAME check_U_V_slices_fail_diff
          COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
          -s ${nc_file} -c "horiz_winds(:,:,1,:)=T_mid(:,:,:)"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (check_U_V_slices_fail_diff PROPERTIES
          WILL_FAIL TRUE
          FIXTURES_REQUIRED shoc_generate_output_nc_files)

# Wrong layout
add_test (NAME check_U_V_slices_fail_layout
          COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
          -s ${nc_file} -c "horiz_winds(:,:,1,:)=U(:,1,:)"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (check_U_V_slices_fail_layout PROPERTIES
          WILL_FAIL TRUE
          FIXTURES_REQUIRED shoc_generate_output_nc_files)

# Missing variable(s)
add_test (NAME check_U_V_slices_fail_missing
          COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
          -s ${nc_file} -c "T=T"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (check_U_V_slices_fail_missing PROPERTIES
          WILL_FAIL TRUE
          FIXTURES_REQUIRED shoc_generate_output_nc_files)
