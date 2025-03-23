include (ScreamUtils)

# Find PyTorch and set up dependencies
# Use a configurable path for PyTorch with a default value
if(NOT DEFINED TORCH_CMAKE_PATH)
  # Default path can be overridden by setting TORCH_CMAKE_PATH
  set(TORCH_CMAKE_PATH "/global/u1/m/mahf708/.micromamba/envs/fme/lib/python3.10/site-packages/torch/share/cmake/Torch" 
      CACHE PATH "Path to PyTorch CMake files")
endif()
set(CMAKE_PREFIX_PATH ${TORCH_CMAKE_PATH})
find_package(Torch REQUIRED)
if(NOT Torch_FOUND)
  message(FATAL_ERROR "PyTorch not found. Please set TORCH_CMAKE_PATH to the directory containing TorchConfig.cmake")
endif()

# Include the main SCREAM libraries
include(EkatSetCompilerFlags)

# Create a unit test with all required dependencies
set(TEST_SRCS
  ${SCREAM_SRC_DIR}/share/util/eamxx_ad_test.cpp
)

add_executable(ace_standalone ${TEST_SRCS})
target_include_directories(ace_standalone PUBLIC ${SCREAM_SRC_DIR})
target_include_directories(ace_standalone SYSTEM PUBLIC ${TORCH_INCLUDE_DIRS})

target_link_libraries(ace_standalone PUBLIC
  ace
  ${TORCH_LIBRARIES}
)

# Set labels and properties
set_target_properties(ace_standalone PROPERTIES
  LABELS "ace;physics"
)

# Set AD configurable options
set (NUM_STEPS 1)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
