include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

# Enable Sandy Bridge arch in Kokkos
option(Kokkos_ARCH_SNB "" ON)

set(CMAKE_CXX_FLAGS "-w -cxxlib=/usr/tce/packages/gcc/gcc-8.3.1/rh" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/gcc/gcc-8.3.1/rh/lib/gcc/x86_64-redhat-linux/8/ -mkl" CACHE STRING "" FORCE)

set(SCREAM_INPUT_ROOT "/usr/gdata/e3sm/ccsm3data/inputdata/" CACHE STRING "")
