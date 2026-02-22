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

# Test for StubBackend
add_executable(test_inference_stub_backend test_inference_stub_backend.cpp)
target_link_libraries(test_inference_stub_backend PRIVATE emulator_common)
target_include_directories(test_inference_stub_backend PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME inference_stub_backend_tests COMMAND test_inference_stub_backend)

# Test for data_view
add_executable(test_data_view test_data_view.cpp)
target_link_libraries(test_data_view PRIVATE emulator_common)
target_include_directories(test_data_view PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME data_view_tests COMMAND test_data_view)

# Test for tensor_adapter
add_executable(test_tensor_adapter test_tensor_adapter.cpp)
target_link_libraries(test_tensor_adapter PRIVATE emulator_common)
target_include_directories(test_tensor_adapter PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME tensor_adapter_tests COMMAND test_tensor_adapter)

# Test for io_adapter
add_executable(test_io_adapter test_io_adapter.cpp)
target_link_libraries(test_io_adapter PRIVATE emulator_common)
target_include_directories(test_io_adapter PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME io_adapter_tests COMMAND test_io_adapter)

# Test for data_source and MctDataSource
add_executable(test_data_source test_data_source.cpp)
target_link_libraries(test_data_source PRIVATE emulator_common)
target_include_directories(test_data_source PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME data_source_tests COMMAND test_data_source)

# Test for data_pipeline (integration)
add_executable(test_data_pipeline test_data_pipeline.cpp)
target_link_libraries(test_data_pipeline PRIVATE emulator_common)
target_include_directories(test_data_pipeline PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME data_pipeline_tests COMMAND test_data_pipeline)
