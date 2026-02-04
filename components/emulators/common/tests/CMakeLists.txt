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
