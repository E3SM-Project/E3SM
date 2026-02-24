cmake_minimum_required(VERSION 3.16)

set(EMULATOR_DRIVER_SOURCES
  src/emulator_factory.cpp
)

add_library(emulator_driver ${EMULATOR_DRIVER_SOURCES})

# This target needs to see the base emulator headers + component headers
target_include_directories(emulator_driver PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../common/src
  ${CMAKE_CURRENT_SOURCE_DIR}/../eatm/src
)

# Link to the common framework + components
target_link_libraries(emulator_driver PUBLIC
  emulator_common
  eatm
)

target_compile_features(emulator_driver PUBLIC cxx_std_17)
add_subdirectory(tests)
