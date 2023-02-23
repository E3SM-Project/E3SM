# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
set(SCREAM_MACHINE "mappy" CACHE STRING "")

include (${EKAT_MACH_FILES_PATH}/mappy.cmake)

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-n" CACHE STRING "")
set(SCREAM_MPI_EXTRA_ARGS "--cpu_bind=threads" CACHE STRING "")
set(SCREAM_MPI_THREAD_FLAG "-c" CACHE STRING "")
