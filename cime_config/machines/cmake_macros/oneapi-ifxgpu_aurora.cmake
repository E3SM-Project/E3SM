
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -fsycl-device-code-split=per_kernel -fsycl-max-parallel-link-jobs=16 -Wl,--no-relax -mcmodel=large -flink-huge-device-code")
if (compile_threaded)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fiopenmp -fopenmp-targets=spir64")
endif()

string(APPEND KOKKOS_OPTIONS " -DCMAKE_CXX_STANDARD=17 -DKokkos_ENABLE_SERIAL=On -DKokkos_ARCH_INTEL_PVC=On -DKokkos_ENABLE_SYCL=On -DKokkos_ENABLE_EXPLICIT_INSTANTIATION=Off -DCMAKE_POSITION_INDEPENDENT_CODE=ON")
string(APPEND SYCL_FLAGS " -\-intel -fsycl -fsycl-targets=spir64_gen -mlong-double-64 -mcmodel=large ")
string(APPEND OMEGA_SYCL_EXE_LINKER_FLAGS " -Xsycl-target-backend \"-device pvc\" ")

# Let's start with the best case: using device buffers in MPI calls by default.
# This is paired with MPIR_CVAR_ENABLE_GPU=1 in config_machines.xml. If this
# ends up causing instability, we can switch to OFF.
set(SCREAM_MPI_ON_DEVICE ON CACHE STRING "")
