include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# Load all kokkos settings from Ekat's mach file
set(Kokkos_ENABLE_SERIAL TRUE CACHE BOOL "")

# Enable Broadwell arch in Kokkos
option(Kokkos_ARCH_BDW "" ON)

set(CMAKE_CXX_FLAGS "-w -cxxlib=/usr/tce/packages/gcc/gcc-8.3.1/rh" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/gcc/gcc-8.3.1/rh/lib/gcc/x86_64-redhat-linux/8/ -mkl" CACHE STRING "" FORCE)

set(SCREAM_INPUT_ROOT "/usr/gdata/climdat/ccsm3data/inputdata" CACHE STRING "")
