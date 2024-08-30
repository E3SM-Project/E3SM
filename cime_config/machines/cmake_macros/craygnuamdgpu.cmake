set(MPICC "cc")
set(MPICXX "hipcc") # Needs MPICH_CXX to use hipcc
set(MPIFC "ftn") # Linker needs to be the Cray wrapper ftn, not mpif90
set(SCC "cc")
set(SCXX "hipcc")
set(SFC "ftn")

string(APPEND CPPDEFS " -DLINUX -DFORTRANUNDERSCORE -DNO_R16 -DCPRGNU -DSCREAM_SYSTEM_WORKAROUND_P3_PART2")
if (COMP_NAME STREQUAL gptl)
    string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CMAKE_Fortran_FLAGS " -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -fallow-argument-mismatch")

string(APPEND CMAKE_C_FLAGS_DEBUG " -O0 -g -Wall -fbacktrace -fcheck=bounds -ffpe-trap=invalid,zero,overflow")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -Wall -fbacktrace -fcheck=bounds -ffpe-trap=zero,overflow")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g -Wall -fbacktrace")

string(APPEND CMAKE_C_FLAGS_RELEASE " -g -O2")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -g -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -g -O2")

if (COMP_NAME STREQUAL csm_share)
  string(APPEND CMAKE_C_FLAGS " -std=c99")
endif()
string(APPEND CMAKE_Fortran_FORMAT_FIXED_FLAG " -ffixed-form")
string(APPEND CMAKE_Fortran_FORMAT_FREE_FLAG " -ffree-form")

set(E3SM_LINK_WITH_FORTRAN "TRUE")
string(APPEND CMAKE_CXX_FLAGS " -I$ENV{MPICH_DIR}/include")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L$ENV{ROCM_PATH}/lib -lamdhip64")

if (compile_threaded)
  string(APPEND CMAKE_C_FLAGS " -fopenmp")
  string(APPEND CMAKE_Fortran_FLAGS " -fopenmp")
  string(APPEND CMAKE_CXX_FLAGS " -fopenmp=libgomp")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fopenmp")
endif()

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On -DKokkos_ENABLE_OPENMP=Off")

set(USE_HIP "TRUE")
string(APPEND CMAKE_HIP_FLAGS "$ENV{CXXFLAGS} --offload-arch=gfx90a -munsafe-fp-atomics")
