#==============================================================================
# Emulator Test Support Library
#
# Provides:
# - Custom Catch2 test session with MPI support
# - Test data generators
# - Common test utilities
#==============================================================================

add_library(emulator_test_support
  emulator_test_session.cpp
)

target_include_directories(emulator_test_support PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Include Catch2 headers
if(DEFINED EMULATOR_CATCH2_INCLUDE_DIR)
  target_include_directories(emulator_test_support PUBLIC 
    ${EMULATOR_CATCH2_INCLUDE_DIR}
  )
endif()

# Link against emulator_common and MPI
target_link_libraries(emulator_test_support PUBLIC
  emulator_common
  MPI::MPI_CXX
)

target_compile_features(emulator_test_support PUBLIC cxx_std_17)

# Define CATCH_CONFIG_RUNNER for custom main (PRIVATE - only for this library)
target_compile_definitions(emulator_test_support PRIVATE
  CATCH_CONFIG_RUNNER
)

#------------------------------------------------------------------------------
# Test Data Directory (in build directory)
#------------------------------------------------------------------------------
set(EMULATOR_TEST_DATA_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/test_data")
file(MAKE_DIRECTORY ${EMULATOR_TEST_DATA_BUILD_DIR})

# Define test data directory path (pointing to build dir, not source)
target_compile_definitions(emulator_test_support PUBLIC
  EMULATOR_TEST_DATA_DIR="${EMULATOR_TEST_DATA_BUILD_DIR}"
)

#------------------------------------------------------------------------------
# Generate Dummy Model for LibTorch Tests
#------------------------------------------------------------------------------
if(EMULATOR_HAS_LIBTORCH)
  find_package(Python3 COMPONENTS Interpreter)
  
  if(Python3_FOUND)
    # Generate dummy model at build time
    set(DUMMY_MODEL_PATH "${EMULATOR_TEST_DATA_BUILD_DIR}/dummy_model.pt")
    set(CREATE_MODEL_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/data/create_dummy_model.py")
    
    add_custom_command(
      OUTPUT ${DUMMY_MODEL_PATH}
      COMMAND ${Python3_EXECUTABLE} ${CREATE_MODEL_SCRIPT} ${DUMMY_MODEL_PATH}
      DEPENDS ${CREATE_MODEL_SCRIPT}
      COMMENT "Generating dummy TorchScript model for tests"
      VERBATIM
    )
    
    add_custom_target(generate_test_model ALL
      DEPENDS ${DUMMY_MODEL_PATH}
    )
    
    message(STATUS "  LibTorch tests will generate dummy model at: ${DUMMY_MODEL_PATH}")
  else()
    message(WARNING "Python3 not found - LibTorch tests may fail (cannot generate dummy model)")
  endif()
endif()

message(STATUS "  Test support library configured")
message(STATUS "  Test data dir: ${EMULATOR_TEST_DATA_BUILD_DIR}")
