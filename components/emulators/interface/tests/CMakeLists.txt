if(EMULATOR_ENABLE_LIBTORCH)
  add_executable(test_harness test_harness.cpp)
  target_link_libraries(test_harness PRIVATE emulator_common emulator_inference Catch2 emulator_interface)
  target_compile_definitions(test_harness PRIVATE
      TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
  add_test(NAME test_harness COMMAND test_harness)
endif()
