if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_SLASHPROC")
endif()
string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2")
set(PIO_FILESYSTEM_HINTS "gpfs")
string(APPEND SLIBS " /lcrc/group/e3sm/soft/improv/netlib-lapack/3.12.0/gcc-12.3.0/liblapack.a /lcrc/group/e3sm/soft/improv/netlib-lapack/3.12.0/gcc-12.3.0/libblas.a")

