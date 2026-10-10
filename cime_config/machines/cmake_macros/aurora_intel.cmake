
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
if (compile_threaded)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fiopenmp -fopenmp-targets=spir64")
endif()
set(Kokkos_ENABLE_SYCL FALSE CACHE BOOL "")

# EAMxx ignores generic CMAKE_CXX_FLAGS, includes CMAKE_CXX_FLAGS_[RELEASE,DEBUG]
string(APPEND CMAKE_CXX_FLAGS         " -fp-model=consistent")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -fp-model=consistent")
string(APPEND CMAKE_CXX_FLAGS_DEBUG   " -fp-model=consistent")

string(APPEND CMAKE_Fortran_FLAGS " -fp-model=consistent")
