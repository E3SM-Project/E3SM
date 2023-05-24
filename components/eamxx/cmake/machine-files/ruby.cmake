include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# Load all kokkos settings from Ekat's mach file
set(Kokkos_ENABLE_SERIAL TRUE CACHE BOOL "")

# Enable Broadwell arch in Kokkos
option(Kokkos_ARCH_BDW "" ON)
# Load Broadwell flags and openmp backend for kokkos
include (${EKAT_MACH_FILES_PATH}/kokkos/intel-bdw.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/ -qmkl" CACHE STRING "" FORCE)
set(SCREAM_INPUT_ROOT "/usr/gdata/climdat/ccsm3data/inputdata" CACHE STRING "")
