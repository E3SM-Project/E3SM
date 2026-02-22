# Common emulator library

add_library(emulator_common
  emulator.cpp
  inference/stub_inference_backend.cpp
  inference/create_inference_backend.cpp
  data/data_view.cpp
  data/raw_tensor_adapter.cpp
  data/mct_data_source.cpp
  data/data_pipeline.cpp
)

target_compile_features(emulator_common PUBLIC cxx_std_17)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/inference
  ${CMAKE_CURRENT_SOURCE_DIR}/data
)
