enable_testing()
# Common tests

# emulator_add_test(<target> <source> <test name> [LIBS <libs>...] [WITH_DATA_DIR])
#   Builds a Catch2 test executable and registers it with ctest.
function(emulator_add_test target source test_name)
  cmake_parse_arguments(ARG "WITH_DATA_DIR" "" "LIBS" ${ARGN})
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE ${ARG_LIBS} Catch2)
  if(ARG_WITH_DATA_DIR)
    target_compile_definitions(${target} PRIVATE
      TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/input-data")
  endif()
  add_test(NAME ${test_name} COMMAND ${target})
endfunction()

emulator_add_test(test_initialization test_emulator.cpp emulator_initialization
  LIBS emulator_common WITH_DATA_DIR)
emulator_add_test(test_inference_config test_inference_config.cpp inference_config_tests
  LIBS emulator_inference emulator_common)
emulator_add_test(test_inference_stub_backend test_inference_stub_backend.cpp
  inference_stub_backend_tests LIBS emulator_inference emulator_common)
emulator_add_test(test_inference_tensor test_inference_tensor.cpp inference_tensor_tests
  LIBS emulator_inference)
emulator_add_test(test_emulator_coupling test_emulator_coupling.cpp emulator_coupling
  LIBS emulator_common WITH_DATA_DIR)
emulator_add_test(test_emulator_config test_emulator_config.cpp test_emulator_config
  LIBS emulator_common WITH_DATA_DIR)
emulator_add_test(test_fpe_guard test_fpe_guard.cpp fpe_guard_tests
  LIBS emulator_inference emulator_common)

# Test for LibTorchBackend
#
# TorchScript files can only be written from Python, so the small modules the
# test loads are generated at build time. Runs on the CPU; a GPU is used only
# if one is present.
if(EMULATOR_ENABLE_LIBTORCH)
  execute_process(COMMAND ${EMULATOR_TORCH_PYTHON} -c "import torch"
    RESULT_VARIABLE EMULATOR_TORCH_PYTHON_FAILED OUTPUT_QUIET ERROR_QUIET)

  if(EMULATOR_TORCH_PYTHON_FAILED)
    message(WARNING "Skipping the LibTorch backend test: "
      "'${EMULATOR_TORCH_PYTHON}' cannot import torch. "
      "Set -DEMULATOR_TORCH_PYTHON=/path/to/python to enable it.")
  else()
    set(EMULATOR_TORCH_FIXTURE_DIR ${CMAKE_CURRENT_BINARY_DIR}/torch_fixtures)
    set(EMULATOR_TORCH_FIXTURES
      ${EMULATOR_TORCH_FIXTURE_DIR}/affine.pt
      ${EMULATOR_TORCH_FIXTURE_DIR}/tuple.pt
      ${EMULATOR_TORCH_FIXTURE_DIR}/list.pt
      ${EMULATOR_TORCH_FIXTURE_DIR}/channels.pt)
    set(EMULATOR_TORCH_FIXTURE_SCRIPT
      ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/make_libtorch_fixture.py)

    add_custom_command(
      OUTPUT ${EMULATOR_TORCH_FIXTURES}
      COMMAND ${EMULATOR_TORCH_PYTHON} ${EMULATOR_TORCH_FIXTURE_SCRIPT}
              ${EMULATOR_TORCH_FIXTURE_DIR}
      DEPENDS ${EMULATOR_TORCH_FIXTURE_SCRIPT}
      COMMENT "Writing the LibTorch test fixtures")
    add_custom_target(libtorch_test_fixtures DEPENDS ${EMULATOR_TORCH_FIXTURES})

    emulator_add_test(test_inference_libtorch_backend test_inference_libtorch_backend.cpp
      inference_libtorch_backend_tests LIBS emulator_inference)
    add_dependencies(test_inference_libtorch_backend libtorch_test_fixtures)
    target_compile_definitions(test_inference_libtorch_backend PRIVATE
      EMULATOR_TEST_TORCH_FIXTURE_DIR="${EMULATOR_TORCH_FIXTURE_DIR}")
  endif()
endif()

# Test for PythonBackend
#
# The embedded interpreter needs numpy at run time. If it also has torch,
# the test runs a torch.nn module as well.
if(EMULATOR_ENABLE_PYTHON)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy"
    RESULT_VARIABLE EMULATOR_NUMPY_MISSING OUTPUT_QUIET ERROR_QUIET)

  if(EMULATOR_NUMPY_MISSING)
    message(WARNING "Skipping the Python backend test: "
      "'${Python3_EXECUTABLE}' cannot import numpy.")
  else()
    emulator_add_test(test_inference_python_backend test_inference_python_backend.cpp
      inference_python_backend_tests LIBS emulator_inference)
    execute_process(COMMAND ${Python3_EXECUTABLE} -c
        "import site; print(':'.join(site.getsitepackages()))"
      OUTPUT_VARIABLE EMULATOR_PYTHON_SITE_PACKAGES
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    target_compile_definitions(test_inference_python_backend PRIVATE
      EMULATOR_TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures"
      EMULATOR_TEST_SITE_PACKAGES="${EMULATOR_PYTHON_SITE_PACKAGES}")

    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import torch"
      RESULT_VARIABLE EMULATOR_PYTHON_TORCH_MISSING OUTPUT_QUIET ERROR_QUIET)
    if(NOT EMULATOR_PYTHON_TORCH_MISSING)
      target_compile_definitions(test_inference_python_backend PRIVATE
        EMULATOR_TEST_HAVE_TORCH)
    endif()

    # A relocatable Python does not find its standard library when embedded
    execute_process(COMMAND ${Python3_EXECUTABLE} -c
        "import sys; print(sys.base_prefix)"
      OUTPUT_VARIABLE EMULATOR_PYTHON_HOME OUTPUT_STRIP_TRAILING_WHITESPACE)
    set_tests_properties(inference_python_backend_tests PROPERTIES
      ENVIRONMENT "PYTHONHOME=${EMULATOR_PYTHON_HOME}")
  endif()
endif()
