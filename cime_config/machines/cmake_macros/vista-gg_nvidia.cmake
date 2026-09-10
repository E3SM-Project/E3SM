string(APPEND CPPDEFS " -DLINUX")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DBIT64 -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY -DHAVE_MPI")
endif()
string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -g")

set(HOMME_QUAD_PREC FALSE CACHE BOOL "") # nvidia does not seem to support QUAD

if (compile_threaded)
  string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_OPENMP=Off") # work-around for nvidia as kokkos is not passing "-mp" for threaded build
endif()

set(MPICC "mpicc")
set(MPICXX "mpicxx")
set(MPIFC "mpif90")
set(SCC "gcc")
set(SCXX "g++")
set(SFC "gfortran")
