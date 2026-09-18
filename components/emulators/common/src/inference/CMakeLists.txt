# Inference library
#
# Depends on nothing else in the emulator framework, so that anything that
# needs to run a network (a component, a single process) can link it alone.

option(EMULATOR_ENABLE_LIBTORCH
  "Build the LibTorch inference backend (needs a LibTorch install)" OFF)

set(EMULATOR_INFERENCE_SOURCES
  tensor.cpp
  fpe_guard.cpp
  inference_config.cpp
  inference_backend.cpp
  stub_inference_backend.cpp
  create_inference_backend.cpp
)

if(EMULATOR_ENABLE_LIBTORCH)
  # Point CMAKE_PREFIX_PATH (or Torch_DIR) at the LibTorch install
  find_package(Torch REQUIRED)
  list(APPEND EMULATOR_INFERENCE_SOURCES libtorch_inference_backend.cpp)
endif()

add_library(emulator_inference ${EMULATOR_INFERENCE_SOURCES})

target_compile_features(emulator_inference PUBLIC cxx_std_17)

target_include_directories(emulator_inference PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if(EMULATOR_ENABLE_LIBTORCH)
  # Recent LibTorch headers need C++20; only the backend's source sees them
  target_compile_features(emulator_inference PRIVATE cxx_std_20)
  target_compile_definitions(emulator_inference PUBLIC EMULATOR_ENABLE_LIBTORCH)
  # PUBLIC: TORCH_LIBRARIES carries the C++ ABI flag LibTorch was built with
  target_link_libraries(emulator_inference PUBLIC "${TORCH_LIBRARIES}")
endif()
