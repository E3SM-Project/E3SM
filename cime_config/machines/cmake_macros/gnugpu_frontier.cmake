set(MPICC "cc")
set(MPICXX "hipcc")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "hipcc")
set(SFC "ftn")

string(APPEND CPPDEFS " -DLINUX")
if (COMP_NAME STREQUAL gptl)
    string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND FFLAGS " -Wno-implicit-interface")

if (NOT DEBUG)
  string(APPEND CFLAGS   " -O2")
  string(APPEND CXXFLAGS " -O2")
  string(APPEND FFLAGS   " -O2")
endif()
string(APPEND CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")

string(APPEND CXXFLAGS " -I$ENV{MPICH_DIR}/include --offload-arch=gfx90a")
string(APPEND SLIBS    " -Wl,--copy-dt-needed-entries -L/opt/cray/pe/gcc-libs -lgfortran -L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa ")

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On")

set(USE_HIP "TRUE")
string(APPEND HIP_FLAGS "${CXXFLAGS} -munsafe-fp-atomics -x hip")
