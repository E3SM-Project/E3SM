set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")

string(APPEND CPPDEFS " -DLINUX")
if (COMP_NAME STREQUAL gptl)
    string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/opt/cray/pe/gcc/12.2.0/snos/lib64 -lgfortran -lstdc++")

# to support Fortran specific compiler intrinsic functions
set(E3SM_LINK_WITH_FORTRAN "TRUE")
