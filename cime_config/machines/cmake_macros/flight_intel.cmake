set(ALBANY_PATH "/projects/ccsm/AlbanyTrilinos_20190904/albany-build/install")
set(CLDERA_PATH "/projects/cldera/cldera-tools/install-master/intel/")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_SLASHPROC")
endif()
set(ESMF_LIBDIR "/projects/ccsm/esmf-6.3.0rp1/lib/libO/Linux.intel.64.openmpi.default")
set(PIO_FILESYSTEM_HINTS "lustre")

# Hack similar to pm-cpu_intel.cmake. Reset the following flags and load
# with oneapi supported options.
set(CMAKE_CXX_FLAGS " ")
set(CMAKE_Fortran_FLAGS_DEBUG " ")

if (compile_threaded)
  string(APPEND CMAKE_CXX_FLAGS " -qopenmp")
endif()
string(APPEND CMAKE_CXX_FLAGS " -fp-model precise")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -O0 -g ")

if (COMP_NAME STREQUAL gcam)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,--no-relax")
  string(APPEND CMAKE_Fortran_FLAGS " -mcmodel=medium")
  string(APPEND CMAKE_C_FLAGS " -mcmodel=medium")
  string(APPEND CMAKE_CXX_FLAGS " -DNDEBUG")
endif()
