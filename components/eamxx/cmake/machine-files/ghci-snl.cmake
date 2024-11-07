# Set Fortran flags
set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch" CACHE STRING "Fortran compiler flags" FORCE)

# Set the path to SCREAM input data
set(SCREAM_INPUT_ROOT /projects/e3sm/inputdata CACHE PATH "Path to SCREAM input data" FORCE)

# Set the path to BLAS/LAPACK libraries
set(BLAS_LIBRARIES "/spack-installs/netlib-lapack/3.11.0/gcc/12.3.0/base/65x6uge/lib64/libblas.so" CACHE STRING "Path to BLAS library" FORCE)
set(LAPACK_LIBRARIES "/spack-installs/netlib-lapack/3.11.0/gcc/12.3.0/base/65x6uge/lib64/liblapack.so" CACHE STRING "Path to LAPACK library" FORCE)

# Let's catch usage of code deprecated in Kokkos 4
option (Kokkos_ENABLE_DEPRECATED_CODE_4 "" OFF)

# We need to manage resources to spread across available cores/gpus
option (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES "" ON)

# Needed by EkatCreateUnitTest
set (EKAT_MPIRUN_EXE "mpirun" CACHE STRING "")
set (EKAT_MPI_NP_FLAG "-n" CACHE STRING "")
