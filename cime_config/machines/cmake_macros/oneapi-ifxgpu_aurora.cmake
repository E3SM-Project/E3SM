
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
if (compile_threaded)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fiopenmp -fopenmp-targets=spir64")
endif()
string(APPEND SYCL_FLAGS " -\-intel -fsycl -fsycl-targets=spir64_gen -mlong-double-64 ")
string(APPEND OMEGA_SYCL_EXE_LINKER_FLAGS " -Xsycl-target-backend \"-device 12.60.7\" ")
string(APPEND CMAKE_CXX_FLAGS " -Xclang -fsycl-allow-virtual-functions")
