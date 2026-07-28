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
target_link_libraries(test_inference_config PRIVATE emulator_common)
target_include_directories(test_inference_config PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME inference_config_tests COMMAND test_inference_config)

# Test for Tensor / TensorMap
add_executable(test_inference_tensor test_inference_tensor.cpp)
target_link_libraries(test_inference_tensor PRIVATE emulator_common)
target_include_directories(test_inference_tensor PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME inference_tensor_tests COMMAND test_inference_tensor)

# Test for StubBackend
add_executable(test_inference_stub_backend test_inference_stub_backend.cpp)
target_link_libraries(test_inference_stub_backend PRIVATE emulator_common)
target_include_directories(test_inference_stub_backend PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME inference_stub_backend_tests COMMAND test_inference_stub_backend)

# Test for InferenceContext
add_executable(test_inference_context test_inference_context.cpp)
target_link_libraries(test_inference_context PRIVATE emulator_common)
target_include_directories(test_inference_context PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME inference_context_tests COMMAND test_inference_context)

# Test for the embedded-Python backend
if(EMULATOR_ENABLE_PYTHON)
  add_executable(test_inference_python_backend test_inference_python_backend.cpp)
  target_link_libraries(test_inference_python_backend PRIVATE emulator_common)
  target_include_directories(test_inference_python_backend PRIVATE ${CATCH2_INCLUDE_DIR})
  target_compile_definitions(test_inference_python_backend PRIVATE
    EMULATOR_TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
  add_test(NAME inference_python_backend_tests COMMAND test_inference_python_backend)
endif()
