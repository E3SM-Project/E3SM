# Standalone unit test of histComputeMod (no CIME, no pFUnit).
#
#   cmake -S . -B build [-DCMAKE_Fortran_COMPILER=nvfortran] [-DOPENMP=ON] [-DACC=MULTICORE|GPU]
#   cmake --build build && ctest --test-dir build --output-on-failure
#   build/test_hist_compute bench     # timing: old algorithm vs histComputeMod
cmake_minimum_required(VERSION 3.16)
project(hist_compute_test Fortran)

option(OPENMP "Thread histComputeMod with OpenMP" OFF)
set(ACC "OFF" CACHE STRING "OpenACC: OFF, MULTICORE or GPU (nvfortran)")

set(E3SM_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../../../../..)
set(ELM_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/../..)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -fbacktrace -Wall -Wno-compare-reals -Wno-target-lifetime")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "^(NVHPC|PGI)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mpreprocess")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Mbounds")
endif()

if (OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

add_executable(test_hist_compute
  ${E3SM_ROOT}/share/util/shr_kind_mod.F90
  ${ELM_MAIN}/histComputeMod.F90
  fake_subgrid_mod.F90
  hist_reference_mod.F90
  test_hist_compute.F90)
if (OPENMP)
  target_link_libraries(test_hist_compute PRIVATE OpenMP::OpenMP_Fortran)
endif()
if (ACC STREQUAL "MULTICORE")
  target_compile_options(test_hist_compute PRIVATE -acc=multicore -Minfo=accel)
  target_link_options(test_hist_compute PRIVATE -acc=multicore)
elseif (ACC STREQUAL "GPU")
  target_compile_options(test_hist_compute PRIVATE -acc=gpu -Minfo=accel)
  target_link_options(test_hist_compute PRIVATE -acc=gpu)
endif()

enable_testing()
add_test(NAME hist_compute COMMAND test_hist_compute)
