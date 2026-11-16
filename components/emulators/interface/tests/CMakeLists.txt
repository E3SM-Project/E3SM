if(EMULATOR_ENABLE_LIBTORCH)
  add_executable(test_harness test_harness.cpp)
  target_link_libraries(test_harness PRIVATE emulator_common emulator_inference Catch2 emulator_interface)
  target_compile_definitions(test_harness PRIVATE
    TEST_DATA_DIR="${CMAKE_CURRENT_BINARY_DIR}/fixtures")
  add_test(NAME test_harness COMMAND test_harness)
  add_dependencies(test_harness emulator_interface_test_fixtures)

  # copy test fixtures into our binary fixtures directory
  set(EMULATOR_INTERFACE_FIXTURES
    ${CMAKE_CURRENT_BINARY_DIR}/fixtures/land_emulator.pt
  )
    
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fixtures DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_command(
    OUTPUT ${EMULATOR_INTERFACE_FIXTURES}
    COMMAND ${CMAKE_COMMAND} -E copy ${EMULATOR_TEST_FIXTURE_DIR}/land_emulator.pt ${CMAKE_CURRENT_BINARY_DIR}/fixtures
    COMMENT "Generating emulator interface test fixtures..."
  )
  add_custom_target(emulator_interface_test_fixtures DEPENDS ${EMULATOR_INTERFACE_FIXTURES})
  add_dependencies(emulator_interface_test_fixtures emulator_coupling_test_fixtures)
endif()
