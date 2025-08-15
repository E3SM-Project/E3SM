set(MPICC "cc")
set(MPICXX "mpicxx")
set(MPIFC "ftn")
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

string(APPEND CMAKE_CXX_FLAGS " --offload-arch=gfx90a")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L$ENV{ROCM_PATH}/lib -lamdhip64")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/opt/cray/pe/gcc/12.2.0/snos/lib64 -lgfortran -lstdc++")
#
#string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/opt/rocm-5.4.0/lib -lhsa-runtime64 -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa ")

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On -DKokkos_ENABLE_OPENMP=Off")

set(USE_HIP "TRUE")
string(APPEND CMAKE_HIP_FLAGS "${CXXFLAGS} -munsafe-fp-atomics -x hip")


set(E3SM_LINK_WITH_FORTRAN "TRUE")
