include (ScreamUtils)

# Set common AD configurable options
set (NUM_STEPS 1)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Create executable for CldFraction standalone
CreateADUnitTestExec(cld_fraction_standalone
  LIBS cld_fraction)

#####################################
#        CldFraction (CPP)          #
#####################################

# Configure yaml files to run directory
set (POSTFIX cpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_cpp.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_cpp.yaml)

CreateUnitTestFromExec(cld_fraction_standalone_cpp cld_fraction_standalone
  EXE_ARGS "--args -ifile=input_cpp.yaml"
  LABELS cld_fraction physics
  FIXTURES_SETUP cldfrac_cpp)

if (EAMXX_ENABLE_PYTHON)
  #####################################
  #        CldFraction (PY)           #
  #####################################

  include (BuildCprnc)
  BuildCprnc()

  # Configure yaml files to run directory
  set (POSTFIX py)
  set (PY_MODULE_NAME "cld_fraction_numpy")
  set (PY_MODULE_PATH ${SCREAM_SOURCE_DIR}/src/physics/cld_fraction)
  set (PY_BACKEND "host")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_py.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_py.yaml)

  # Test the process with python impl
  CreateUnitTestFromExec(cld_fraction_standalone_py cld_fraction_standalone
    EXE_ARGS "--args -ifile=input_py.yaml"
    LABELS cld_fraction physics
    FIXTURES_SETUP cldfrac_py)

  # Compare output of py and cpp tests
  set (SRC_FILE "cldfrac_standalone_output_cpp.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
  set (TGT_FILE "cldfrac_standalone_output_py.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
  set (TEST_NAME cldfrac_standalone_cpp_vs_py)
  add_test (NAME ${TEST_NAME}
            COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TEST_NAME} PROPERTIES
        LABELS "cldfrac;infrastructure"
        FIXTURES_REQUIRED "cldfrac_py;cldfrac_cpp")

  if (Kokkos_ENABLE_CUDA)
    # Also run with cupy instead of numpy
    set (PY_MODULE_NAME "cld_fraction_cupy")
    set (PY_BACKEND "device")
    set (POSTFIX "cupy")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/input_cupy.yaml)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/output_cupy.yaml)

    # Test the process with cupy impl
    CreateUnitTestFromExec(cld_fraction_standalone_cupy cld_fraction_standalone
      EXE_ARGS "--args -ifile=input_cupy.yaml"
      LABELS cld_fraction physics
      FIXTURES_SETUP cldfrac_cupy)

    # Compare output of cupy and cpp tests
    set (SRC_FILE "cldfrac_standalone_output_cpp.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
    set (TGT_FILE "cldfrac_standalone_output_cupy.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
    set (TEST_NAME cldfrac_standalone_cpp_vs_cupy)
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES
          LABELS "cldfrac;infrastructure"
          FIXTURES_REQUIRED "cldfrac_cupy;cldfrac_cpp")
  endif()

  if (NOT EAMXX_ENABLE_GPU)
    #####################################
    #         CldFracNet (PY)           #
    #####################################

    # Create executable that runs cld_frac_net emulator
    CreateADUnitTestExec(cld_frac_net_standalone
      LIBS cld_frac_net)

    # Test the process with python ml emulator
    set (PY_MODULE_NAME "cld_frac_net")
    set (PY_MODULE_PATH ${SCREAM_BASE_DIR}/src/physics/cld_fraction/cld_frac_net)
    set (POSTFIX py)
    set (EMULATOR pytorch)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_ml.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/input_ml_${POSTFIX}.yaml)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_ml.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/output_ml_${POSTFIX}.yaml)

    CreateUnitTestFromExec(cld_frac_net_${POSTFIX} cld_frac_net_standalone
      EXE_ARGS "--args -ifile=input_ml_${POSTFIX}.yaml"
      LABELS cld_fraction physics
      FIXTURES_SETUP cld_frac_net_${POSTFIX})

    # So that we know, below, that we can run cpp-vs-py test for cld_frac_net output
    set (HAS_CLD_FRAC_NET_PY TRUE)
  endif()
endif()


if (NOT EAMXX_ENABLE_GPU)
  #####################################
  #         CldFracNet (CPP)          #
  #####################################

  # Create test that runs lapis-generated emulator
  unset (PY_MODULE_NAME)
  set (POSTFIX cpp)
  set (EMULATOR lapis)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_ml.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_ml_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_ml.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_ml_${POSTFIX}.yaml)

  CreateUnitTestFromExec(cld_frac_net_${POSTFIX} cld_fraction_standalone
    EXE_ARGS "--args -ifile=input_ml_${POSTFIX}.yaml"
    LABELS cld_fraction
    FIXTURES_SETUP cld_frac_net_${POSTFIX})

  if (HAS_CLD_FRAC_NET_PY)
    # Compare output of py and cpp tests of cld_frac_net
    set (SRC_FILE "cld_frac_net_standalone_output_cpp.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
    set (TGT_FILE "cld_frac_net_standalone_output_py.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
    set (TEST_NAME cld_frac_net_standalone_cpp_vs_py)
    add_test (NAME ${TEST_NAME}
      COMMAND ${SCREAM_BASE_DIR}/scripts/compare-nc-files
                -s ${SRC_FILE} -t ${TGT_FILE} --tol 1e-6
                -c cldfrac_ice=cldfrac_ice cldfrac_tot=cldfrac_tot
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES
          LABELS "cldfrac;infrastructure"
          FIXTURES_REQUIRED "cld_frac_net_py;cld_frac_net_cpp")
  endif()
endif()
