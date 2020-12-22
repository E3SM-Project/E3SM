cmake_minimum_required(VERSION 3.0)
project(standalone)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()
enable_testing()

add_subdirectory(${YAKL_HOME} ./yakl)
include_directories(${YAKL_HOME})

file(GLOB CPP_SRC     *.cpp)

add_executable(test_fft_cpp.exe ${CPP_SRC})

target_link_libraries(test_fft_cpp.exe yakl )

set_property(TARGET test_fft_cpp.exe APPEND PROPERTY COMPILE_FLAGS ${DEFS2D} )

