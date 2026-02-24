# Common source library
cmake_minimum_required(VERSION 3.16)

set(EMULATOR_COMMON_CXX_SOURCES
  emulator.cpp
  emulator_c_api.cpp
  inference/stub_inference_backend.cpp
  inference/create_inference_backend.cpp
)

set(EMULATOR_COMMON_F90_SOURCES
  emulator_f2c_api.F90
  emulator_f_api.F90
  emulator_handle.F90
)

add_library(emulator_common
  ${EMULATOR_COMMON_CXX_SOURCES}
  ${EMULATOR_COMMON_F90_SOURCES}
)

target_compile_features(emulator_common PUBLIC cxx_std_17)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/inference
)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
)

