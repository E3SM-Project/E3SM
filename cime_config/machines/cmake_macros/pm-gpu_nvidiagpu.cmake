string(APPEND CONFIG_ARGS " --host=cray")
set(USE_CUDA "TRUE")
set(HOMME_QUAD_PREC FALSE CACHE BOOL "") # nvidia does not support QUAD (real*16)
string(APPEND CPPDEFS " -DGPU")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CMAKE_Fortran_FLAGS " -tp=zen3")
string(APPEND CMAKE_C_FLAGS       " -tp=zen3")
string(APPEND CMAKE_CXX_FLAGS     " -tp=zen3")
# Disable OpenACC (use CUDA via Kokkos instead).
string(APPEND CMAKE_C_FLAGS " -noacc")
string(APPEND CMAKE_Fortran_FLAGS " -noacc")
string(APPEND CMAKE_CXX_FLAGS " -noacc")
if (COMP_NAME STREQUAL cice)
  # nvfortran -Mbounds produces false-positive subscript-out-of-range errors
  # for module-level fixed-size arrays (e.g. TLON) accessed inside OpenMP
  # parallel regions in CICE.  Disable bounds checking for CICE only.
  string(APPEND CMAKE_Fortran_FLAGS " -Mnobounds")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -Mnobounds")
endif()
string(APPEND CMAKE_CUDA_FLAGS " -ccbin CC -O2 -arch sm_80 --use_fast_math")
string(APPEND KOKKOS_OPTIONS " -DKokkos_ARCH_AMPERE80=On -DKokkos_ENABLE_CUDA=On -DKokkos_ENABLE_CUDA_LAMBDA=On -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=Off -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF -DCMAKE_CXX_FLAGS=-noacc -DCMAKE_C_FLAGS=-noacc")
set(CMAKE_CUDA_ARCHITECTURES "80")
# NVHPC embeds __acc_compiled in every object file unless compiled with -noacc.
# EKAT objects (e.g. libekat_yamlparser) may reference this symbol, so we must
# satisfy it at link time via libacchost.
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_cuda -noacc -lacchost")
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")
