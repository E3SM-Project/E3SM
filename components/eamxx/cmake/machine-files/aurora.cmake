include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake)

set(EKAT_MPIRUN_EXE "mpiexec" CACHE STRING "" FORCE)
set(EKAT_MPI_NP_FLAG "-np" CACHE STRING "" FORCE)
set(EKAT_MPI_EXTRA_ARGS "--label --cpu-bind depth -envall" CACHE STRING "")
set(EKAT_MPI_THREAD_FLAG "-d" CACHE STRING "")

set(NETCDF_C_PATH "$ENV{NETCDF_C_PATH}")
set(NETCDF_FORTRAN_PATH "$ENV{NETCDF_FORTRAN_PATH}")
set(PNETCDF_PATH "$ENV{PNETCDF_PATH}")
