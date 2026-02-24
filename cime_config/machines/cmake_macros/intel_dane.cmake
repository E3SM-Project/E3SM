string(APPEND CPPDEFS " -DNO_SHR_VMATH -DCNL")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check all -ftrapuv")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/")
set(KOKKOS_OPTIONS "--with-serial --ldflags='-L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/'")
