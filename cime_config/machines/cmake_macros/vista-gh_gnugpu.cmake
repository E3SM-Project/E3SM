set(USE_CUDA "TRUE")
string(APPEND CPPDEFS " -DGPU")
string(APPEND CPPDEFS " -DLINUX")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DBIT64 -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY -DHAVE_MPI")
  # -DHAVE_NANOTIME -- cant use this as the assembly instructions that wont work on ARM
endif()
string(APPEND CMAKE_CUDA_FLAGS " -ccbin CC -O2 -arch sm_90 --use_fast_math")
string(APPEND KOKKOS_OPTIONS " -DKokkos_ARCH_HOPPER90=On -DKokkos_ENABLE_CUDA=On -DKokkos_ENABLE_CUDA_LAMBDA=On -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=Off -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=Off")
set(CMAKE_CUDA_ARCHITECTURES "90")

string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_EXE_LINKER_FLAGS "  -L/opt/apps/gcc/14.2.0/lib64 -lstdc++") # workaround for pnetcdf thinking it needs gcc14 abi

set(MPICC "mpicc")
set(MPICXX "mpicxx") # Needs MPICH_CXX to use nvcc <env name="MPICH_CXX">$SHELL{which nvcc}</env>
set(MPIFC "mpif90")
set(SCC "gcc")
set(SCXX "nvcc")
set(SFC "gfortran")
