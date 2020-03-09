# Note: CMAKE_CXX_COMPILER needs to be set to the path of nvcc_wrapper
# nvcc_wrapper will choose either the Nvidia Cuda compiler or the OpenMP compiler depending on what's being compiled

# TPL settings
set(NETCDF_DIR $ENV{NETCDF_ROOT} CACHE FILEPATH "")
set(NetCDF_Fortran_PATH /ascldap/users/lbertag/workdir/libs/netcdf/netcdf-f/netcdf-f-install/waterman/gcc CACHE FILEPATH "")
set(HDF5_DIR $ENV{HDF5_ROOT} CACHE FILEPATH "")
set(ZLIB_DIR $ENV{ZLIB_ROOT} CACHE FILEPATH "")
set(CURL_ROOT $ENV{CURL_ROOT} CACHE FILEPATH "")
set(CURL_LIBRARY -L$ENV{CURL_ROOT}/lib -lcurl CACHE LIST "")

# Flag tweaks
set(CMAKE_C_FLAGS "-w" CACHE STRING "")
set(CMAKE_CXX_FLAGS "-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -Wno-unknown-pragmas --fmad=false -G" CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-w" CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS "-ldl" CACHE STRING "")
set(OPT_FLAGS "-O0" CACHE STRING "")
set(DEBUG_FLAGS "-ffp-contract=off -g"CACHE STRING "")

# Homme settings
set(AVX_VERSION 0 CACHE STRING "")
set(HOMMEXX_VECTOR_SIZE 1 CACHE STRING "")
set(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
set(USE_QUEUING FALSE CACHE BOOL "")
set(USE_TRILINOS OFF CACHE BOOL "")
set(WITH_PNETCDF FALSE CACHE FILEPATH "")
set(USE_NUM_PROCS 1 CACHE STRING "Num mpiprocs to use")
set(HAVE_EXTRAE TRUE CACHE BOOL "")
set(Extrae_LIBRARY "-L${NETCDF_DIR}/lib -L${NetCDF_Fortran_PATH}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -ldl -lz" CACHE STRING "")
set(ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set(ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")
set(HOMMEXX_EXEC_SPACE "CUDA" CACHE STRING "")

# Kokkos settings
set(ENABLE_OPENMP FALSE CACHE BOOL "")
set(KOKKOS_ENABLE_DEBUG FALSE CACHE BOOL "")
set(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION FALSE CACHE BOOL "")
set(KOKKOS_ENABLE_CUDA TRUE CACHE BOOL "")
set(KOKKOS_ENABLE_CUDA_LAMBDA TRUE CACHE BOOL "")
set(KOKKOS_ARCH Volta70 CACHE STRING "")
set(KOKKOS_ENABLE_DEPRECATED_CODE FALSE CACHE BOOL "")
set(KOKKOS_ENABLE_EXPLICIT_INSTANTIATION FALSE CACHE BOOL "")

# Compilers
set(CMAKE_C_COMPILER mpicc CACHE STRING "")
set(CMAKE_Fortran_COMPILER mpif90 CACHE STRING "")
set(CMAKE_CXX_COMPILER mpicxx CACHE STRING "")
