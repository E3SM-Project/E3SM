
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
if (compile_threaded)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fiopenmp -fopenmp-targets=spir64")
endif()
set(Kokkos_ENABLE_SYCL FALSE CACHE BOOL "")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.3")
  # ifx: clear out -freciprocal-math/unsafe-fp-math from the default fast model; source first clears them
  string(PREPEND CMAKE_Fortran_FLAGS " -fp-model=source")
endif()

