include (ScreamUtils)

# P3's warm-rain emulator is a python module, run through the python inference
# backend of components/emulators
if (NOT EAMXX_ENABLE_INFERENCE)
  return()
endif()

# The SDM-trained emulator needs torch in the python the backend embeds
execute_process(COMMAND ${EAMXX_INFERENCE_PYTHON} -c "import numpy, torch"
  RESULT_VARIABLE NO_TORCH OUTPUT_QUIET ERROR_QUIET)
if (NO_TORCH)
  # (an error, not a skip: the test must not go missing unnoticed)
  message(FATAL_ERROR "Test p3_warm_rain_emulator needs numpy and torch, but '${EAMXX_INFERENCE_PYTHON}' "
    "cannot import them. Set -DPython3_EXECUTABLE=/path/to/python.")
endif()

set (TEST_BASE_NAME p3_warm_rain_emulator)
set (ATM_TIME_STEP 100)
set (NUM_STEPS 6)
set (RUN_T0 2021-10-12-45000)

CreateADUnitTestExec(${TEST_BASE_NAME}
  LIBS p3)

GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})

# The models: synthetic ones (synthetic_emulators.py), and the SDM-trained emulator
set (EMULATOR_DIR ${SCREAM_SRC_DIR}/physics/p3/warm_rain_emulator)
set (PYTHON_PATH "${CMAKE_CURRENT_SOURCE_DIR}:${EMULATOR_DIR}")

# Creates test p3_warm_rain_emulator_${CASE}, on RANKS MPI ranks
macro (CreateCase CASE USE_EMULATOR MODULE MODEL FALLBACK_SCALE RANKS)
  # Macro arguments are not variables: make them visible to configure_file
  set (CASE ${CASE})
  set (USE_EMULATOR ${USE_EMULATOR})
  set (MODULE ${MODULE})
  set (MODEL ${MODEL})
  set (FALLBACK_SCALE ${FALLBACK_SCALE})
  # An optional last argument: whether the emulator replaces P3's rates
  set (REPLACES_P3 false)
  if (${ARGC} GREATER 6)
    set (REPLACES_P3 ${ARGV6})
  endif()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml ${CMAKE_CURRENT_BINARY_DIR}/input_${CASE}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml ${CMAKE_CURRENT_BINARY_DIR}/output_${CASE}.yaml)
  CreateUnitTestFromExec(${TEST_BASE_NAME}_${CASE} ${TEST_BASE_NAME}
    EXE_ARGS "--args -ifile=input_${CASE}.yaml"
    MPI_RANKS ${RANKS}
    LABELS p3 physics
    FIXTURES_SETUP_INDIVIDUAL ${TEST_BASE_NAME}_${CASE})
endmacro()

# Compares the outputs of two cases, within the optional tolerance (default 0: BFB)
macro (CompareCases NAME CASE1 RANKS1 CASE2 RANKS2)
  set (TOL 0)
  if (${ARGC} GREATER 5)
    set (TOL ${ARGV5})
  endif()
  add_test (NAME ${TEST_BASE_NAME}_${NAME}
    COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
    -s ${TEST_BASE_NAME}_${CASE1}.INSTANT.nsteps_x1.np${RANKS1}.${RUN_T0}.nc
    -t ${TEST_BASE_NAME}_${CASE2}.INSTANT.nsteps_x1.np${RANKS2}.${RUN_T0}.nc --tol ${TOL}
    -c T_mid=T_mid qv=qv qc=qc qr=qr qi=qi qm=qm nc=nc nr=nr ni=ni bm=bm
       precip_liq_surf_mass=precip_liq_surf_mass precip_ice_surf_mass=precip_ice_surf_mass
       qc2qr_autoconv=qc2qr_autoconv qc2qr_accret=qc2qr_accret
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TEST_BASE_NAME}_${NAME} PROPERTIES
    LABELS "p3;physics"
    FIXTURES_REQUIRED "${TEST_BASE_NAME}_${CASE1}_np${RANKS1}_omp1;${TEST_BASE_NAME}_${CASE2}_np${RANKS2}_omp1")
endmacro()

# P3 without the emulator
CreateCase(stock false none none 1.0 1)

# An emulator that never applies: the whole path runs (inputs, model, merge),
# and P3 must be BFB with P3 without it
CreateCase(nowhere true synthetic_emulators nowhere 1.0 1)
CompareCases(nowhere_vs_stock nowhere 1 stock 1)

# What the emulator returns is what P3 uses, and what P3 gives it is its state:
# - an emulator of zero rates that applies everywhere: no autoconversion, no accretion;
# - an emulator that converts all the cloud water it is given to rain: no cloud water left.
macro (CreateCheckedCase CASE MODEL CHECK)
  CreateCase(${CASE} true synthetic_emulators ${MODEL} 1.0 1 ${ARGN})
  add_test (NAME ${TEST_BASE_NAME}_${CASE}_check
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/check_synthetic.py ${CHECK}
            ${TEST_BASE_NAME}_${CASE}.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
            ${TEST_BASE_NAME}_stock.INSTANT.nsteps_x1.np1.${RUN_T0}.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TEST_BASE_NAME}_${CASE}_check PROPERTIES
    LABELS "p3;physics"
    FIXTURES_REQUIRED "${TEST_BASE_NAME}_${CASE}_np1_omp1;${TEST_BASE_NAME}_stock_np1_omp1")
endmacro()
CreateCheckedCase(no_rain no_rain no_rain)
CreateCheckedCase(drain drain drain)

# An emulator that replaces P3's rates is used in every cell, whatever its masks
# say: the one that never applies (and returns large rates) leaves no cloud water
CreateCheckedCase(replace nowhere drain true)

# The SDM-trained emulator, with the autoconversion fallback its authors suggest.
# torch does not always compute a grid point in exactly the same way in batches of
# different sizes, so runs on different numbers of ranks agree to round-off of the
# (single precision) network only.
set (RANKS 1)
if (TEST_RANK_END GREATER 1)
  list (APPEND RANKS ${TEST_RANK_END})
endif()
CreateCase(sdm true sdm_warm_rain ${EMULATOR_DIR}/sdm_warm_rain.pt 0.36 "${RANKS}")
if (TEST_RANK_END GREATER 1)
  CompareCases(sdm_np1_vs_np${TEST_RANK_END} sdm 1 sdm ${TEST_RANK_END} 1e-4)
endif()
