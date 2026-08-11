string(APPEND CONFIG_ARGS " --host=cray")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()

set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "icx")
set(SCXX "icpx")
set(SFC "ifx")

# CPU-only Intel oneAPI build on GPU nodes (no CUDA/Kokkos GPU settings).
# Same approach as pm-cpu_intel.cmake: reset CXX flags to drop -fp-model=source
# (unsupported by icpx), then re-apply the correct flags.
set(CMAKE_CXX_FLAGS " ")
if (compile_threaded)
  string(APPEND CMAKE_CXX_FLAGS " -qopenmp")
endif()
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g")

# Check for Intel LLVM (ifx) version 2025 or newer
if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
        string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check nouninit") # Applying Intel 2025.3 sanitization workaround
    endif()
endif()

string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CXX_FLAGS " -fp-model=precise")
string(APPEND CMAKE_CXX_FLAGS " -fp-model=consistent")
string(APPEND CMAKE_Fortran_FLAGS " -fp-model=consistent -fimf-use-svml")
string(APPEND CMAKE_Fortran_FLAGS " -DHAVE_ERF_INTRINSICS")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -g -traceback")
