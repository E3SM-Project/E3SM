# Common tests
# Test for emulator_registry
add_executable(test_emulator_registry test_emulator_registry.cpp)
target_link_libraries(test_emulator_registry PRIVATE emulator_common Catch2)
add_test(NAME emulator_registry_tests COMMAND test_emulator_registry)

# Test for emulator
add_executable(test_emulator test_emulator.cpp)
target_link_libraries(test_emulator PRIVATE emulator_common Catch2)
add_test(NAME emulator_tests COMMAND test_emulator)

# Test for InferenceConfig and utilities
add_executable(test_inference_config test_inference_config.cpp)
target_link_libraries(test_inference_config PRIVATE emulator_common Catch2)
add_test(NAME inference_config_tests COMMAND test_inference_config)

# Test for StubBackend
add_executable(test_inference_stub_backend test_inference_stub_backend.cpp)
target_link_libraries(test_inference_stub_backend PRIVATE emulator_common Catch2)
add_test(NAME inference_stub_backend_tests COMMAND test_inference_stub_backend)


add_executable(test_ekat_yaml
    test_yaml_fields.cpp
)

target_link_libraries(test_ekat_yaml
    PRIVATE
        ekat::YamlParser
        Catch2
)
target_compile_definitions(test_ekat_yaml
    PRIVATE
        TEST_FIELDS_YAML="${CMAKE_CURRENT_SOURCE_DIR}/test_fields.yaml"
)

add_test(
    NAME test_ekat_yaml
    COMMAND test_ekat_yaml
)


