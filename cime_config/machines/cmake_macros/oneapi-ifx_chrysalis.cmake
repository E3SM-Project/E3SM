set(MPIFC "mpif90")
set(MPICC "mpicc")
set(MPICXX "mpicxx")
string(APPEND SLIBS " -fiopenmp -fopenmp-targets=spir64")
