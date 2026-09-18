# Inference library
#
# Depends on nothing else in the emulator framework, so that anything that
# needs to run a network (a component, a single process) can link it alone.

set(EMULATOR_INFERENCE_SOURCES
  tensor.cpp
  inference_config.cpp
  inference_backend.cpp
  stub_inference_backend.cpp
  create_inference_backend.cpp
)

add_library(emulator_inference ${EMULATOR_INFERENCE_SOURCES})

target_compile_features(emulator_inference PUBLIC cxx_std_17)

target_include_directories(emulator_inference PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
