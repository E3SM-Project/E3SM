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

set(HOMME_QUAD_PREC FALSE CACHE BOOL "") # nvidia does not support QUAD (real*16)

string(APPEND CMAKE_C_FLAGS          " -tp=grace")
string(APPEND CMAKE_CXX_FLAGS        " -tp=grace")
string(APPEND CMAKE_Fortran_FLAGS    " -tp=grace")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -tp=grace -nonvmalloc")
