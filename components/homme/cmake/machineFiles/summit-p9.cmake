#interactive job
#bsub -W 2:00 -nnodes 1 -P cli115 -Is /bin/bash


#cmake -C ~/acme-fork-lb/components/homme/cmake/machineFiles/summit.cmake -DHOMMEXX_MPI_ON_DEVICE=FALSE ~/acme-fork-lb/components/homme/

#SET (HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")

#modules, note no cuda, otherwise CUDA_BUILD=true
#module load cmake/3.6.1  gcc/6.4.0 netlib-lapack/3.6.1
#module load netcdf/4.6.1 netcdf-fortran/4.4.4
#module load hdf5/1.10.3


SET (NETCDF_DIR $ENV{OLCF_NETCDF_FORTRAN_ROOT} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{OLCF_HDF5_ROOT} CACHE FILEPATH "")

SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(USE_TRILINOS OFF CACHE BOOL "")
SET(E3SM_KOKKOS_PATH "/ccs/home/onguba/kokkos/build-serial-p9-nodebug" CACHE STRING "")

SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "mpicc" CACHE STRING "")
#SET(CMAKE_CXX_COMPILER "/ccs/home/onguba/kokkos/bin/nvcc_wrapper" CACHE STRING "")

SET (ADD_Fortran_FLAGS "-fopenmp -O3" CACHE STRING "")
SET (ADD_C_FLAGS       "-fopenmp -O3 --std=c++11" CACHE STRING "")
SET (ADD_CXX_FLAGS     "-fopenmp -O3 --std=c++11" CACHE STRING "")

#it is really only for stdc++, not for anything else
SET (HAVE_EXTRAE TRUE CACHE BOOL "")
SET (Extrae_LIBRARY "-L${KOKKOS_PATH}/lib -lkokkos -ldl -lstdc++" CACHE STRING "")

set (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
#set to false to test
set (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

set (HOMME_TESTING_PROFILE "dev" CACHE STRING "")

set (USE_NUM_PROCS 4 CACHE STRING "")

set (OPT_FLAGS " -mcpu=power9 -mtune=power9 -DNDEBUG " CACHE STRING "")
SET (USE_MPI_OPTIONS "--bind-to core" CACHE FILEPATH "")
