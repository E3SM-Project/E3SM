include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

#serial is needed, but maybe it is always on?
include (${EKAT_MACH_FILES_PATH}/kokkos/mi250.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/hip.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

SET(MPICH_DIR "/opt/cray/pe/mpich/8.1.16/ofi/crayclang/10.0" CACHE STRING "")

set(CMAKE_CXX_FLAGS "--amdgpu-target=gfx90a -fno-gpu-rdc  -I${MPICH_DIR}/include" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "-I${MPICH_DIR}/include" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "-I${MPICH_DIR}/include" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.16/gtl/lib -lmpi_gtl_hsa" CACHE STRING "" FORCE)



