
#### DO NOT FORGET
#export LD_LIBRARY_PATH=/home/onguba/kokkos-build/bld-sycl-may04/lib64/:$LD_LIBRARY_PATH


#interactive job


#cmake -C ~/acme-fork-lb/components/homme/cmake/machineFiles/summit.cmake -DHOMMEXX_MPI_ON_DEVICE=FALSE ~/acme-fork-lb/components/homme/

#cmake -C ~/acme-MASTER-GB/components/homme/cmake/machineFiles/crusher-gpumpi.cmake -DE3SM_KOKKOS_PATH=/ccs/home/onguba/kokkos-crusher-june2022/bld-hipcc ~/acme-MASTER-GB/components/homme/

#SET (HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")
SET (HOMMEXX_CUDA_MAX_WARP_PER_TEAM "16" CACHE STRING  "")

#SET (NETCDF_DIR $ENV{OLCF_NETCDF_FORTRAN_ROOT} CACHE FILEPATH "")
#SET (HDF5_DIR $ENV{OLCF_HDF5_ROOT} CACHE FILEPATH "")

SET(BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")

#SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

#temp hack
SET(HOMME_USE_KOKKOS TRUE CACHE BOOL "")

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")

set(KOKKOS_HOME "/home/onguba/kokkos-build/bld-sycl-may04" CACHE STRING "")
set(E3SM_KOKKOS_PATH $ENV{KOKKOS_HOME} CACHE STRING "")

#SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

SET(USE_TRILINOS OFF CACHE BOOL "")

#CUDA_BUILD is set in SetCompilersFlags, after findPackage(Cuda)
#i haven't extend it to hip, set it here instead
SET(SYCL_BUILD TRUE CACHE BOOL "")
SET(HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")

#uncomment this if using internal kokkos build
#SET(Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
####SET(CMAKE_CXX_STANDARD "14" CACHE STRING "")
#SET(Kokkos_ENABLE_DEBUG OFF CACHE BOOL "")
#SET(Kokkos_ARCH_VEGA90A ON CACHE BOOL "")
#SET(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
#SET(Kokkos_ENABLE_HIP ON CACHE BOOL "")
####SET(Kokkos_ENABLE_CUDA_LAMBDA OFF CACHE BOOL "")
#SET(Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")

#set(E3SM_KOKKOS_PATH "/home/onguba/kokkos-build/bld-sycl" CACHE STRING "")
#bld3 is with unnamed-lambdas flag
#bld4 is PVC=off, GEN=on

#set(KOKKOS_HOME "/home/onguba/kokkos-build/bld-sycl4" CACHE STRING "")

SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "")

#set(extra "-Wno-deprecated-declarations -fno-finite-math-only -O3 -DNDEBUG -fsycl -fno-sycl-id-queries-fit-in-int -fsycl-dead-args-optimization -fsycl-unnamed-lambda -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.60.7\" ")
set(extra "-Wno-deprecated-declarations -fno-finite-math-only -O3 -DNDEBUG -fsycl -fno-sycl-id-queries-fit-in-int -fsycl-dead-args-optimization -fsycl-unnamed-lambda ")
#set(extra "-fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.60.7\" ")
#set(extra "-fsycl-link-huge-device-code -fPIC")
#set(extra "-fsycl-link-huge-device-code")

set(extra2 "-L/soft/restricted/CNDA/updates/2022.10.15.001/oneapi/compiler/trunk-20221014/compiler/linux/lib")

SET(ADD_Fortran_FLAGS "-fc=ifx -O3 -DNDEBUG -DCPRINTEL -g" CACHE STRING "")
SET(ADD_C_FLAGS "-cc=icx -O3 -DNDEBUG -g" CACHE STRING "")
#fopenmp flag is for compiler error saying kokkos install has openmp space
#SET(ADD_CXX_FLAGS "-cxx=icpx -std=c++17 -O3 -qopenmp -fsycl ${extra} -DNDEBUG -I$ENV{KOKKOS_HOME}/include" CACHE STRING "")
SET(ADD_CXX_FLAGS "-cxx=icpx -g -std=c++17 -O3 -fsycl ${extra} -DNDEBUG -I${KOKKOS_HOME}/include" CACHE STRING "")
#SET(ADD_LINKER_FLAGS "-L$ENV{KOKKOS_HOME}/lib64 ${extra2} -lsycl" CACHE STRING "")
SET(ADD_LINKER_FLAGS "-L${KOKKOS_HOME}/lib64 ${extra2} -Wno-deprecated-declarations -fno-finite-math-only -O3 -DNDEBUG -DKOKKOS_DEPENDENCE -fsycl -fno-sycl-id-queries-fit-in-int -fsycl-dead-args-optimization " CACHE STRING "")


set (ENABLE_OPENMP OFF CACHE BOOL "")
set (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

set (HOMME_TESTING_PROFILE "dev" CACHE STRING "")

set (USE_NUM_PROCS 4 CACHE STRING "")

#set (OPT_FLAGS "-mcpu=power9 -mtune=power9" CACHE STRING "")
SET (USE_MPI_OPTIONS "--bind-to core" CACHE FILEPATH "")
