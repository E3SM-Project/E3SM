cmake_minimum_required(VERSION 3.0)
project(standalone)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
enable_testing()

file(GLOB FORTRAN_SRC ../fortran/*.F90 ../fortran/SGS_TKE/*.F90 ../fortran/MICRO_SAM1MOM/*.F90 ../fortran/ADV_MPDATA/*.F90 ../fortran/*.c             )
file(GLOB CPP_SRC     ../*.F90     ../SGS_TKE/*.F90     ../MICRO_SAM1MOM/*.F90     ../ADV_MPDATA/*.F90     ../*.c     ../*.cpp)
file(GLOB CUDA_SRC    ../*.cpp)

add_subdirectory(${YAKL_HOME} ./yakl  )
include_directories(${YAKL_HOME})

add_subdirectory(fortran2d)
add_subdirectory(fortran3d)
add_subdirectory(cpp2d)
add_subdirectory(cpp3d)


