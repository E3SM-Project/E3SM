# Common tests

# Test for emulator_context
add_executable(test_emulator_context test_emulator_context.cpp)
target_link_libraries(test_emulator_context PRIVATE emulator_common)
target_include_directories(test_emulator_context PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME emulator_context_tests COMMAND test_emulator_context)

# Test for emulator_comp
add_executable(test_emulator_comp test_emulator_comp.cpp)
target_link_libraries(test_emulator_comp PRIVATE emulator_common)
target_include_directories(test_emulator_comp PRIVATE ${CATCH2_INCLUDE_DIR})
add_test(NAME emulator_comp_tests COMMAND test_emulator_comp)
