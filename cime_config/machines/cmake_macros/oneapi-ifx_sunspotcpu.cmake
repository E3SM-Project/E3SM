
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
if (compile_threaded)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fiopenmp -fopenmp-targets=spir64")
endif()

string(APPEND KOKKOS_OPTIONS " -DCMAKE_CXX_STANDARD=17 -DKokkos_ENABLE_SERIAL=On -DKokkos_ENABLE_EXPLICIT_INSTANTIATION=Off")

#set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "")










