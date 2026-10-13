
string(APPEND CMAKE_EXE_LINKER_FLAGS " -fsycl-device-code-split=per_kernel -fsycl-max-parallel-link-jobs=16 -Wl,--no-relax")
if (compile_threaded)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fiopenmp -fopenmp-targets=spir64")
endif()

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_SERIAL=On -DKokkos_ARCH_INTEL_PVC=On -DKokkos_ENABLE_SYCL=On -DKokkos_ENABLE_EXPLICIT_INSTANTIATION=Off -DCMAKE_POSITION_INDEPENDENT_CODE=ON")
string(APPEND SYCL_FLAGS " -fsycl -fsycl-targets=spir64_gen -mlong-double-64 ")
string(APPEND OMEGA_SYCL_EXE_LINKER_FLAGS " -Xsycl-target-backend \"-device pvc\" ")

# Let's start with the best case: using device buffers in MPI calls by default.
# This is paired with MPIR_CVAR_ENABLE_GPU=1 in config_machines.xml. If this
# ends up causing instability, we can switch to OFF.
set(SCREAM_MPI_ON_DEVICE ON CACHE STRING "")

set(USE_SYCL "TRUE")

# Override the value TRUE set by intelgpu.cmake (via intel.cmake)
set(E3SM_LINK_WITH_FORTRAN "FALSE")

# EAMxx ignores generic CMAKE_CXX_FLAGS, includes CMAKE_CXX_FLAGS_[RELEASE,DEBUG]
string(APPEND CMAKE_CXX_FLAGS         " -fp-model=consistent")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -fp-model=consistent")
string(APPEND CMAKE_CXX_FLAGS_DEBUG   " -fp-model=consistent")

string(APPEND CMAKE_Fortran_FLAGS " -fp-model=consistent")

# 'just' -g may lead to linker internal errors and/or huge builds out of quotas
string(APPEND CMAKE_C_FLAGS_DEBUG   " -fno-system-debug")
string(APPEND CMAKE_CXX_FLAGS_DEBUG   " -fno-system-debug")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " --offload-compress")

