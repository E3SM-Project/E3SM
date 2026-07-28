# Common source library
cmake_minimum_required(VERSION 3.16)

add_subdirectory(coupler_api)

set(EMULATOR_COMMON_CXX_SOURCES
  emulator.cpp
  emulator_c_api.cpp
  inference/tensor.cpp
  inference/inference_config.cpp
  inference/inference_context.cpp
  inference/inference_backend.cpp
  inference/stub_inference_backend.cpp
  inference/create_inference_backend.cpp
)

if(EMULATOR_ENABLE_PYTHON)
  list(APPEND EMULATOR_COMMON_CXX_SOURCES
    inference/python_interpreter.cpp
    inference/python_inference_backend.cpp
  )
endif()

set(EMULATOR_COMMON_F90_SOURCES
  emulator_f2c_api.F90
  emulator_handle.F90
)

add_library(emulator_common
  ${EMULATOR_COMMON_CXX_SOURCES}
  ${EMULATOR_COMMON_F90_SOURCES}
)

target_compile_features(emulator_common PUBLIC cxx_std_20)

# Place generated Fortran .mod files in a predictable directory so that
# downstream Fortran targets (e.g. the driver test) can find them.
set_target_properties(emulator_common PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/inference
)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_link_libraries(emulator_common PUBLIC
  E3SM::coupler_api
)

if(EMULATOR_HAVE_MPI)
  target_compile_definitions(emulator_common PUBLIC EMULATOR_HAVE_MPI)
  target_link_libraries(emulator_common PUBLIC MPI::MPI_CXX)
endif()

if(EMULATOR_ENABLE_PYTHON)
  target_compile_definitions(emulator_common PUBLIC EMULATOR_ENABLE_PYTHON)
  target_link_libraries(emulator_common PUBLIC Python3::Python)
  # Where the emulator python package lives, so a run does not have to set
  # PYTHONPATH.  Anything the user names in `inference.python_path` still wins.
  target_compile_definitions(emulator_common PRIVATE
    EMULATOR_PYTHON_PACKAGE_DIR="${EMULATOR_PYTHON_PACKAGE_DIR}")
endif()
