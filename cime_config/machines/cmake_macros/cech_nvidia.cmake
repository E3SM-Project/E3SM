set(MPICC "mpicc")
set(MPICXX "mpicxx")
set(MPIFC "mpifort")
set(SCC "nvc")
set(SCXX "nvc++")
set(SFC "nvfortran")

set(HOMME_QUAD_PREC FALSE CACHE BOOL "") # nvidia does not support QUAD (real*16)

string(APPEND CMAKE_C_FLAGS       " -tp=grace")
string(APPEND CMAKE_CXX_FLAGS     " -tp=grace")
string(APPEND CMAKE_Fortran_FLAGS " -tp=grace")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -tp=grace")
