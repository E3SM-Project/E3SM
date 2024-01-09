if (compile_threaded)
  #string(APPEND CFLAGS " -fopenmp")
  string(APPEND CMAKE_Fortran_FLAGS " -fopenmp")
  string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fopenmp")
endif()
if (COMP_NAME STREQUAL elm)
  string(APPEND CMAKE_Fortran_FLAGS " -hfp0")
endif()
string(APPEND CMAKE_Fortran_FLAGS " -hipa0 -hzero -hsystem_alloc -f free -N 255 -h byteswapio")

string(APPEND CMAKE_EXE_LINKER_FLAGS " -L$ENV{ROCM_PATH}/lib -lamdhip64 $ENV{OLCF_LIBUNWIND_ROOT}/lib/libunwind.a /sw/frontier/spack-envs/base/opt/cray-sles15-zen3/clang-14.0.0-rocm5.2.0/gperftools-2.10-6g5acp4pcilrl62tddbsbxlut67pp7qn/lib/libtcmalloc.a")

set(PIO_FILESYSTEM_HINTS "gpfs")
