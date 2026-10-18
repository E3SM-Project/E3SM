# Standalone unit test of histComputeMod and subgridAveViewMod (no CIME, no pFUnit).
#
#   cmake -S . -B build [-DCMAKE_Fortran_COMPILER=nvfortran] [-DACC=MULTICORE|GPU] [-DGPU_FLAGS=cc80]
#   cmake --build build && ctest --test-dir build --output-on-failure
#   build/test_hist_compute bench     # timing: original algorithm vs histComputeMod
cmake_minimum_required(VERSION 3.16)
project(hist_compute_test Fortran)

set(ACC "OFF" CACHE STRING "OpenACC: OFF, MULTICORE or GPU (nvfortran)")
set(GPU_FLAGS "cc80" CACHE STRING "nvfortran -gpu= suboptions for ACC=GPU")

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
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
endif()

add_executable(test_hist_compute
  ${E3SM_ROOT}/share/util/shr_kind_mod.F90
  ${ELM_MAIN}/subgridAveViewMod.F90
  ${ELM_MAIN}/histComputeMod.F90
  fake_subgrid_mod.F90
  hist_reference_mod.F90
  test_hist_compute.F90)
if (ACC STREQUAL "MULTICORE")
  target_compile_options(test_hist_compute PRIVATE -acc=multicore -Minfo=accel)
  target_link_options(test_hist_compute PRIVATE -acc=multicore)
elseif (ACC STREQUAL "GPU")
  # The reference runs on the host: without FMA on the device, the results are
  # bit-for-bit (nvfortran 26.5 ignores -gpu=nofma at -O0; -Mnofma works at any -O)
  target_compile_options(test_hist_compute PRIVATE -acc=gpu -gpu=${GPU_FLAGS} -Mnofma -Minfo=accel)
  target_link_options(test_hist_compute PRIVATE -acc=gpu -gpu=${GPU_FLAGS})
endif()

enable_testing()
add_test(NAME hist_compute COMMAND test_hist_compute)
