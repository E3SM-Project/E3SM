string(APPEND CPPDEFS " -DLINUX")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check all -ftrapuv")
set(PIO_FILESYSTEM_HINTS "lustre")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lpmi")
