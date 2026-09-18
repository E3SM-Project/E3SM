# Common tests

# Test for emulator_registry
add_executable(test_emulator_registry test_emulator_registry.cpp)
target_link_libraries(test_emulator_registry PRIVATE emulator_common)
target_include_directories(test_emulator_registry PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME emulator_registry_tests COMMAND test_emulator_registry)

# Test for emulator
add_executable(test_emulator test_emulator.cpp)
target_link_libraries(test_emulator PRIVATE emulator_common)
target_include_directories(test_emulator PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME emulator_tests COMMAND test_emulator)

# Test for InferenceConfig and utilities
add_executable(test_inference_config test_inference_config.cpp)
target_link_libraries(test_inference_config PRIVATE emulator_inference)
target_include_directories(test_inference_config PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME inference_config_tests COMMAND test_inference_config)

# Test for StubBackend
add_executable(test_inference_stub_backend test_inference_stub_backend.cpp)
target_link_libraries(test_inference_stub_backend PRIVATE emulator_inference)
target_include_directories(test_inference_stub_backend PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME inference_stub_backend_tests COMMAND test_inference_stub_backend)

# Test for Tensor and TensorMap
add_executable(test_inference_tensor test_inference_tensor.cpp)
target_link_libraries(test_inference_tensor PRIVATE emulator_inference)
target_include_directories(test_inference_tensor PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME inference_tensor_tests COMMAND test_inference_tensor)

# Test for FpeGuard
add_executable(test_fpe_guard test_fpe_guard.cpp)
target_link_libraries(test_fpe_guard PRIVATE emulator_inference)
target_include_directories(test_fpe_guard PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME fpe_guard_tests COMMAND test_fpe_guard)

# Test for LibTorchBackend
#
# TorchScript files can only be written from Python, so the small modules the
# test loads are generated at build time. Runs on the CPU; a GPU is used only
# if one is present.
if(EMULATOR_ENABLE_LIBTORCH)
  set(EMULATOR_TORCH_PYTHON "python3" CACHE STRING
    "Python interpreter with torch, used to write the LibTorch test fixtures")
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

    add_executable(test_inference_libtorch_backend test_inference_libtorch_backend.cpp)
    add_dependencies(test_inference_libtorch_backend libtorch_test_fixtures)
    target_link_libraries(test_inference_libtorch_backend PRIVATE emulator_inference)
    target_include_directories(test_inference_libtorch_backend PRIVATE ${CATCH2_INCLUDE_DIR})
    target_compile_definitions(test_inference_libtorch_backend PRIVATE
      EMULATOR_TEST_TORCH_FIXTURE_DIR="${EMULATOR_TORCH_FIXTURE_DIR}")
    add_test(NAME inference_libtorch_backend_tests COMMAND test_inference_libtorch_backend)
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
    add_executable(test_inference_python_backend test_inference_python_backend.cpp)
    target_link_libraries(test_inference_python_backend PRIVATE emulator_inference)
    target_include_directories(test_inference_python_backend PRIVATE ${CATCH2_INCLUDE_DIR})
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
    add_test(NAME inference_python_backend_tests COMMAND test_inference_python_backend)

    # A relocatable Python does not find its standard library when embedded
    execute_process(COMMAND ${Python3_EXECUTABLE} -c
        "import sys; print(sys.base_prefix)"
      OUTPUT_VARIABLE EMULATOR_PYTHON_HOME OUTPUT_STRIP_TRAILING_WHITESPACE)
    set_tests_properties(inference_python_backend_tests PROPERTIES
      ENVIRONMENT "PYTHONHOME=${EMULATOR_PYTHON_HOME}")
  endif()
endif()
