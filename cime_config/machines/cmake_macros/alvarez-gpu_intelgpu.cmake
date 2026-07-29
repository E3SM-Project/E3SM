string(APPEND CONFIG_ARGS " --host=cray")
set(USE_CUDA "TRUE")
string(APPEND CPPDEFS " -DGPU")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")
string(APPEND KOKKOS_OPTIONS " -DKokkos_ARCH_AMPERE80=On -DKokkos_ENABLE_CUDA=On -DKokkos_ENABLE_CUDA_LAMBDA=On -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=Off -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF")
# Workaround for Intel oneAPI 2025.3+ (clang/21) + NVCC 12.9 incompatibility:
# Intel clang/21 headers use builtins (e.g., __builtin_elementwise_popcount) that NVCC's
# internal Clang (~LLVM 14) does not support. This sets CMAKE_CUDA_HOST_COMPILER=g++ for
# direct cmake-managed .cu file compilation.
#
# NOTE: This does NOT fix the EKAT kokkoscore build failure, which occurs when
# kokkos_launch_compiler routes C++ files through nvcc_wrapper with icpx as ccbin.
# That failure is fixed by pre-setting NVCC_WRAPPER_DEFAULT_COMPILER in share/build/buildlib.ekat.
find_program(GCC_CXX_COMPILER "g++")
if (GCC_CXX_COMPILER)
  string(APPEND CMAKE_CUDA_FLAGS " -ccbin ${GCC_CXX_COMPILER}")
  string(APPEND KOKKOS_OPTIONS " -DCMAKE_CUDA_HOST_COMPILER=${GCC_CXX_COMPILER}")
else()
  message(WARNING "g++ not found; Intel oneAPI 2025.3+ CUDA host compiler incompatibility may occur")
endif()
string(APPEND CMAKE_CUDA_FLAGS " -O2 -arch sm_80 --use_fast_math --allow-unsupported-compiler")
set(CMAKE_CUDA_ARCHITECTURES "80")
set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "icx")
set(SCXX "icpx")
set(SFC "ifx")

# The Intel LLVM icpx on this system does not support -fp-model=source (set by intelgpu.cmake).
# Reset CXX flags and replicate the correct set, same approach as alvarez-cpu_intel.cmake.
set(CMAKE_CXX_FLAGS " ")
if (compile_threaded)
  string(APPEND CMAKE_CXX_FLAGS " -qopenmp")
endif()
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CXX_FLAGS " -fp-model=precise")
string(APPEND CMAKE_CXX_FLAGS " -fp-model=consistent")

# Check for Intel LLVM (ifx) version 2025 or newer
if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
        string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check nouninit") # Applying Intel 2025.3 sanitization workaround
    endif()
endif()

string(APPEND CMAKE_Fortran_FLAGS " -fp-model=consistent -fimf-use-svml")
string(APPEND CMAKE_Fortran_FLAGS " -DHAVE_ERF_INTRINSICS")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -g -traceback")
