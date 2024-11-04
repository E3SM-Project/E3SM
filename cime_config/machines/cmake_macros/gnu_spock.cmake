string(APPEND CMAKE_Fortran_FLAGS " -Wno-implicit-interface ")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
set(PIO_FILESYSTEM_HINTS "gpfs")
