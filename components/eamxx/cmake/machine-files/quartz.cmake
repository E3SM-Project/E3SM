# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
set(Kokkos_ENABLE_SERIAL TRUE CACHE BOOL "")

# Enable Broadwell arch in Kokkos
option(Kokkos_ARCH_BDW "" ON)

#if COMPILER is not defined, should be running standalone with quartz-intel or quartz-gcc
if ("${COMPILER}" STREQUAL "intel")
   set(CMAKE_CXX_FLAGS "-w -cxxlib=/usr/tce/packages/gcc/gcc-8.3.1/rh" CACHE STRING "" FORCE)
   set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/gcc/gcc-8.3.1/rh/lib/gcc/x86_64-redhat-linux/8/ -mkl" CACHE STRING "" FORCE)
elseif ("${COMPILER}" STREQUAL "gnu")
   message(WARNING "You are using an unsupported e3sm compiler. For supported quartz compilers run ./${E3SM_ROOT}/cime/scripts/query_config --machines quartz")
   set(CMAKE_CXX_FLAGS "-w" CACHE STRING "" FORCE)
   set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/gcc/gcc-8.3.1/rh/lib/gcc/x86_64-redhat-linux/8/" CACHE STRING "" FORCE)
endif()

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-n" CACHE STRING "")
set(SCREAM_MPI_EXTRA_ARGS "" CACHE STRING "")
set(SCREAM_INPUT_ROOT "/usr/gdata/climdat/ccsm3data/inputdata" CACHE STRING "")
