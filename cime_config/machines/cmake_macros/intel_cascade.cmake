string(APPEND CONFIG_ARGS " --enable-filesystem-hints=lustre")
string(APPEND CPPDEFS " -DLINUX")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check all -ftrapuv")
set(PIO_FILESYSTEM_HINTS "lustre")
