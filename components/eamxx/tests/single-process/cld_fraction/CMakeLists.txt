include (ScreamUtils)

# Set AD configurable options
set (NUM_STEPS 1)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Configure yaml files to run directory
set (POSTFIX cpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_cpp.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_cpp.yaml)

# Test the process
CreateADUnitTest(cld_fraction_standalone
  EXE_ARGS "--args -ifile=input_cpp.yaml"
  LIBS cld_fraction
  LABELS cld_fraction physics
  FIXTURES_SETUP cldfrac_cpp)

if (EAMXX_ENABLE_PYTHON)
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
    # Run an ml emulator for cld-fraction
    set (PY_MODULE_NAME "cld_fraction_ml")
    set (PY_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
    set (POSTFIX pyml)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/input_pyml.yaml)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/output_pyml.yaml)

    # Test the process with python ml emulator
    CreateUnitTestFromExec(cld_fraction_standalone_pyml cld_fraction_standalone
      EXE_ARGS "--args -ifile=input_pyml.yaml"
      LABELS cld_fraction physics
      FIXTURES_SETUP cldfrac_pyml)
  endif()
endif()
