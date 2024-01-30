set(MPICC "cc")
set(MPICXX "hipcc") # Needs MPICH_CXX to use hipcc
set(MPIFC "ftn") # Linker needs to be the Cray wrapper ftn, not mpif90
set(SCC "cc")
set(SCXX "hipcc")
set(SFC "ftn")

string(APPEND CPPDEFS " -DLINUX")
if (COMP_NAME STREQUAL gptl)
    string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CMAKE_Fortran_FLAGS " -Wno-implicit-interface")

string(APPEND CMAKE_C_FLAGS_RELEASE   " -O2")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE   " -O2")
string(APPEND SPIO_CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")

set(E3SM_LINK_WITH_FORTRAN "TRUE")
string(APPEND CMAKE_CXX_FLAGS " -I$ENV{MPICH_DIR}/include --offload-arch=gfx90a")
string(APPEND CMAKE_EXE_LINKER_FLAGS    " -L/opt/cray/pe/gcc-libs -lgfortran -L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa ")

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On -DKokkos_ENABLE_OPENMP=Off")

set(USE_HIP "TRUE")
string(APPEND CMAKE_HIP_FLAGS "${CXXFLAGS} -munsafe-fp-atomics -x hip")
