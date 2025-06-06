# ==============================================================================
# from Luca's gnu_morellino.cmake
# ==============================================================================

set(NetCDF_Fortran_PATH "$ENV{NetCDF_Fortran_ROOT}" CACHE PATH "Path to netcdf Fortran installation")
set(NETCDF_FORTRAN_PATH "$ENV{NetCDF_Fortran_ROOT}" CACHE PATH "")
set(NetCDF_Fortran_ROOT "$ENV{NetCDF_Fortran_ROOT}" CACHE PATH "")
set(NetCDF_Fortran_LIBRARY "$ENV{NetCDF_Fortran_ROOT}/lib/libnetcdff.so" CACHE FILEPATH "")
set(netcdf_fortran_lib "$ENV{NetCDF_Fortran_ROOT}/lib/libnetcdff.so" CACHE FILEPATH "")
set(NetCDF_Fortran_INCLUDE_DIR "$ENV{NetCDF_Fortran_ROOT}/include" CACHE PATH "" FORCE)
set(NetCDF_Fortran_INCLUDE_DIRS "$ENV{NetCDF_Fortran_ROOT}/include" CACHE PATH "" FORCE)
set(NetCDF_C_PATH "$ENV{NetCDF_C_ROOT}" CACHE PATH "")
set(NETCDF_C_PATH "$ENV{NetCDF_C_ROOT}" CACHE PATH "")
set(PnetCDF_C_PATH "$ENV{PnetCDF_C_ROOT}" CACHE PATH "")
set(HDF5_PATH "$ENV{HDF_ROOT}" CACHE PATH "")
# set(BLAS_PATH "$ENV{lapack_path}" CACHE PATH "")
# set(BLAS_LIBRARIES "$ENV{lapack_path}/lib64/libblas.so" CACHE FILEPATH "")
# set(LAPACK_LIBRARIES "$ENV{lapack_path}/lib64/liblapack.so" CACHE FILEPATH "")
# if (MPILIB STREQUAL mpi-serial)
#   string(APPEND FFLAGS " -I/usr/lib64/gfortran/modules")
# endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17 -Wno-error=template-body -Wno-error=expansion-to-defined -Wno-variadic-macros -Wno-varargs"  CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu17"  CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans"  CACHE STRING "" FORCE)
# string(APPEND FFLAGS " -fallow-argument-mismatch")

# Load all kokkos settings from Ekat's mach file
set(EKAT_MACH_FILES_PATH "$ENV{HOME}/Sandia/research/eamxx_project/e3sm2/externals/ekat/cmake/machine-files")
include(${EKAT_MACH_FILES_PATH}/kokkos/generic.cmake)

# set(SCREAM_INPUT_ROOT "$ENV{SCREAM_BASELINES_DIR}" CACHE PATH "")

option(EKAT_SKIP_FIND_YAML_CPP
  "Skip find_package for yaml-cpp, and  build directly from submodule" ON)
option(EKAT_ENABLE_FPE "Turn off FPEs for MacOS build" OFF)

# set(EKAT_MPIRUN_EXE "mpiexec" CACHE STRING "The executable name for mpirun")
# set(EKAT_MPI_EXTRA_ARGS "--bind-to core" CACHE STRING "Options for mpirun")
# set(EKAT_MPI_NP_FLAG "--map-by" CACHE STRING "The mpirun flag for designating the total number of ranks")
# set(EKAT_MPI_THREAD_FLAG "" CACHE STRING "The mpirun flag for designating the number of threads")


# Set testing options
set(SCREAM_TEST_MAX_THREADS 4 CACHE STRING "Upper limit on threads for threaded tests")
set(SCREAM_TEST_THREAD_INC  1 CACHE STRING "Thread count increment for threaded tests")
set(SCREAM_TEST_MAX_RANKS   4 CACHE STRING "Upper limit on ranks for mpi tests")

# Set inputs/baselines paths
# set(SCREAM_INPUT_ROOT $ENV{WORK_DIR}/e3sm/e3sm-data/inputdata CACHE PATH "Path to scream input files")
set(SCREAM_INPUT_ROOT $ENV{SCREAM_INPUT_DIR} CACHE PATH "Path to scream input files")
set(SCREAM_BASELINES_DIR $ENV{SCREAM_BASELINES_DIR}/full_debug CACHE PATH "Master baselines dir")
set(EKAT_ENABLE_TESTS OFF CACHE BOOL "Whether to enable EKAT tests")

# Set TPLs paths
# set(NetCDF_C_PATH $ENV{NETCDF_C_ROOT} CACHE PATH "Path to netcdf C installation")
# set(NetCDF_Fortran_PATH $ENV{NETCDF_F_ROOT} CACHE PATH "Path to netcdf Fortran installation")
# set(BLAS_LIBRARIES   $ENV{BLAS_ROOT}/lib64/libblas.so   CACHE FILEPATH "Path to BLAS libraries")
# set(LAPACK_LIBRARIES $ENV{BLAS_ROOT}/lib64/liblapack.so CACHE FILEPATH "Path to LAPACK libraries")

option (TEST_LAUNCHER_MANAGE_RESOURCES "" ON)
