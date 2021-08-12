cmake_minimum_required(VERSION 3.0)
project(standalone)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${YAKL_ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()
enable_testing()

file(GLOB FORTRAN_SRC ../../sam/*.F90 ../../sam/SGS_TKE/*.F90 ../../sam/MICRO_SAM1MOM/*.F90 ../../sam/ADV_MPDATA/*.F90 ../../sam/*.c )
file(GLOB CPP_SRC     ../*.F90 ../*.cpp)
file(GLOB CUDA_SRC    ../*.cpp)

set(YAKL_BIN ${CMAKE_CURRENT_BINARY_DIR}/yakl)
add_subdirectory(${YAKL_HOME} ./yakl)

add_subdirectory(fortran2d)
add_subdirectory(fortran3d)
add_subdirectory(cpp2d)
add_subdirectory(cpp3d)


