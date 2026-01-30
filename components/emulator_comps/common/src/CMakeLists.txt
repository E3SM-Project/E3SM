# Common source library

add_library(emulator_common
  emulator_comp.cpp
)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
