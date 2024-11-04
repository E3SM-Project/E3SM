include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# Load all kokkos settings from Ekat's mach file
set(Kokkos_ENABLE_SERIAL TRUE CACHE BOOL "")

# Enable Broadwell arch in Kokkos
option(Kokkos_ARCH_BDW "" ON)

#if COMPILER is not defined, should be running standalone with quartz-intel or quartz-gcc
if ("${COMPILER}" STREQUAL "intel")
   set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/mkl/mkl-2022.1.0/lib/intel64/ -qmkl" CACHE STRING "" FORCE)
elseif ("${COMPILER}" STREQUAL "gnu")
   message(WARNING "You are using an unsupported e3sm compiler. For supported quartz compilers run ./${E3SM_ROOT}/cime/scripts/query_config --machines quartz")
   set(CMAKE_CXX_FLAGS "-w" CACHE STRING "" FORCE)
   set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/gcc/gcc-8.3.1/rh/lib/gcc/x86_64-redhat-linux/8/" CACHE STRING "" FORCE)
endif()

set(SCREAM_INPUT_ROOT "/usr/gdata/e3sm/ccsm3data/inputdata" CACHE STRING "")
