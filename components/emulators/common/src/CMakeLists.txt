# Common source library

find_package(MPI REQUIRED)

add_library(emulator_common
  emulator.cpp
)

target_compile_features(emulator_common PUBLIC cxx_std_17)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(emulator_common PUBLIC MPI::MPI_CXX)

