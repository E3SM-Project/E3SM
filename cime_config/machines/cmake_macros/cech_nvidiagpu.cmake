execute_process(COMMAND uname -m
  OUTPUT_VARIABLE CECH_HOST_ARCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if (NOT CECH_HOST_ARCH STREQUAL "aarch64")
  message(FATAL_ERROR
    "On cech, builds must run on aarch64 compute nodes (detected '${CECH_HOST_ARCH}' login node). "
    "Pass '--batch-build' to create_test/case.build (or use batch-create-test.sh), or run inside a compute node allocation."
  )
endif()

set(MPICC "mpicc")
set(MPICXX "mpicxx")
set(MPIFC "mpifort")
set(SCC "nvc")
set(SCXX "nvc++")
set(SFC "nvfortran")

set(USE_CUDA "TRUE")
set(HOMME_QUAD_PREC FALSE CACHE BOOL "") # nvidia does not support QUAD (real*16)
string(APPEND CPPDEFS " -DGPU")
# Revisit once Hommexx launch-bounds tuning (HOMMEXX_CUDA_MIN_BLOCKS) is in upstream source:
# string(APPEND CPPDEFS " -DHOMMEXX_CUDA_MIN_BLOCKS=3")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")

string(APPEND CMAKE_C_FLAGS       " -tp=grace")
string(APPEND CMAKE_CXX_FLAGS     " -tp=grace")
string(APPEND CMAKE_Fortran_FLAGS " -tp=grace")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -tp=grace")

if (MPAS_OPENACC)
  # OpenACC is enabled for Fortran only in MPAS components (mpas-ocean, mpas-seaice).
  # Keep -noacc for C/C++ to prevent NVHPC from embedding __acc_compiled
  # in C/C++ objects, which causes compilation failures when -mp is also used.
  string(APPEND CMAKE_C_FLAGS   " -noacc")
  string(APPEND CMAKE_CXX_FLAGS " -noacc")
  if (COMP_NAME MATCHES "^mpas" OR CORE STREQUAL "ocean" OR CORE STREQUAL "seaice")
    string(APPEND CMAKE_Fortran_FLAGS " -acc -gpu=cc100 -Minfo=accel")
    string(APPEND CPPDEFS " -DMPAS_OPENACC")
  else()
    string(APPEND CMAKE_Fortran_FLAGS " -noacc")
  endif()
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -acc -gpu=cc100 -Minfo=accel")
else()
  # Default Kokkos/CUDA behavior
  string(APPEND CMAKE_C_FLAGS       " -noacc")
  string(APPEND CMAKE_CXX_FLAGS     " -noacc")
  string(APPEND CMAKE_Fortran_FLAGS " -noacc")
  # NVHPC embeds __acc_compiled in every object file unless compiled with -noacc.
  # EKAT objects (e.g. libekat_yamlparser) may reference this symbol, so we must
  # satisfy it at link time via libacchost.
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -noacc -lacchost")
endif()

string(APPEND CMAKE_CUDA_FLAGS " -ccbin mpicxx -O3 -arch sm_100 --use_fast_math")
if (CMAKE_CXX_FLAGS_RELEASE MATCHES "-O2")
  string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
elseif (NOT CMAKE_CXX_FLAGS_RELEASE MATCHES "-O3")
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O3")
endif()
if (DEFINED CMAKE_CUDA_FLAGS_RELEASE AND CMAKE_CUDA_FLAGS_RELEASE MATCHES "-O2")
  string(REPLACE "-O2" "-O3" CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE}")
endif()
string(APPEND KOKKOS_OPTIONS " -DKokkos_ARCH_BLACKWELL100=On -DKokkos_ENABLE_CUDA=On -DKokkos_ENABLE_CUDA_LAMBDA=On -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=Off -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF -DCMAKE_CXX_FLAGS=-noacc -DCMAKE_C_FLAGS=-noacc")
set(KOKKOS_USE_EKAT_MACH_FILE "FALSE")
set(CMAKE_CUDA_ARCHITECTURES "100")
string(REPLACE "-Ktrap=fp" "-Ktrap=divz,ovf" CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
if (COMP_NAME STREQUAL cice)
  string(APPEND CMAKE_Fortran_FLAGS " -Mnobounds")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -Mnobounds")
endif()
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/global/homes/n/ndk/e3sm/3rdparty/compute/hdf5/lib -lhdf5_hl -lhdf5")
