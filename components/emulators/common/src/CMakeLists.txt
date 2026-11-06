# Common source library

set(EMULATOR_COMMON_CXX_SOURCES
  emulator.cpp
  emulator_config.cpp
)


add_library(emulator_common
  ${EMULATOR_COMMON_CXX_SOURCES}
)

target_link_libraries(emulator_common PUBLIC MOAB)
target_link_libraries(emulator_common PRIVATE ekat::YamlParser)
target_compile_features(emulator_common PUBLIC cxx_std_20)

# Place generated Fortran .mod files in a predictable directory so that
# downstream Fortran targets (e.g. the driver test) can find them.
set_target_properties(emulator_common PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_link_libraries(emulator_common PUBLIC emulator_inference)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/inference
)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_link_libraries(emulator_common PUBLIC E3SM::coupler_api)

