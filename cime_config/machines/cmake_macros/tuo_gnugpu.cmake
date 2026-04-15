string(APPEND CONFIG_ARGS " --host=cray")
string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")

set(MPICC "cc")
set(MPICXX "mpicxx") # Needs MPICH_CXX=hipcc
set(MPIFC "ftn") # Linker needs to be the Cray wrapper ftn, not mpif90
set(SCC "cc")
set(SCXX "mpicxx") # kokkos build fail at one point, maybe ok now
set(SFC "ftn")

if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()

string(APPEND CMAKE_Fortran_FLAGS " -Wno-implicit-interface")

string(APPEND SLIBS " -L$ENV{CRAY_NETCDF_HDF5PARALLEL_PREFIX} -L$ENV{CRAY_PARALLEL_NETCDF_PREFIX}/lib -lpnetcdf -lnetcdf -lnetcdff")
string(APPEND SLIBS " -lblas -llapack")
string(APPEND SLIBS " -lxpmem  -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa -Wl,-rpath,$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib") 

string(APPEND SPIO_CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF -DKokkos_ARCH_AMD_GFX942=ON -DKokkos_ARCH_AMD_GFX942_APU=OFF")
set(USE_HIP "TRUE")

string(APPEND CMAKE_HIP_FLAGS "${CXXFLAGS} -O2 -x hip")

string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2")

string(APPEND CXX_LIBS " -lstdc++")

set(E3SM_LINK_WITH_FORTRAN "TRUE")

set(PIO_FILESYSTEM_HINTS "lustre")
set(SCREAM_ENABLE_MAM OFF) # ndk
