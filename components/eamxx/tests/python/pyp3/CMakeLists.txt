include (ScreamUtils)

set (TEST_BASE_NAME p3_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)  # 1h 2.5h 24h
set (RUN_T0 2021-10-12-45000)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)

##########################################
#             Pure CXX run               #
##########################################

set (EXEC_TYPE "cxx")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_cxx.yaml)
CreateADUnitTest(${TEST_BASE_NAME}_cxx
  LABELS p3 physics
  LIBS p3
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}_cxx
)

##########################################
#           Run from python              #
##########################################

set (EXEC_TYPE "py")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_py.yaml)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/p3_standalone_py
                ${CMAKE_CURRENT_BINARY_DIR}/p3_standalone_py
                @ONLY)

# This test will a) run pyp3, and b) compare the generated nc files
CreateUnitTestFromExec(
  ${TEST_BASE_NAME}_py "p3_standalone_py"
  LABELS p3 physics
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}_py
)


CreateRange(MpiRanks ${TEST_RANK_START} ${TEST_RANK_END} 1)
foreach (rank IN LISTS MpiRanks)
  set(suffix "_np${rank}_omp1")
  CompareNCFiles(
    TEST_NAME pyp3_cxx_vs_py_np${rank}
    SRC_FILE p3_standalone_cxx.INSTANT.nsteps_x1.np${rank}.2021-10-12-45000.nc
    TGT_FILE p3_standalone_py.INSTANT.nsteps_x1.np${rank}.2021-10-12-45000.nc
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_cxx${suffix} ${FIXTURES_BASE_NAME}_py${suffix}
    LABELS pyeamxx
  )
endforeach()
