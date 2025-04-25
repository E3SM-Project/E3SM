string(APPEND CONFIG_ARGS " --host=cray")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -g")

# currently, there is known issue with nvidia compiler installation (not seeing all relevant include files)
# and this is temporary work-around github.com/E3SM-Project/E3SM/issues/7003
string(APPEND CMAKE_CXX_FLAGS_RELEASE " --gcc-toolchain=/usr/bin/gcc")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " --gcc-toolchain=/usr/bin/gcc")

if (compile_threaded)
  string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_OPENMP=Off") # work-around for nvidia as kokkos is not passing "-mp" for threaded build
endif()
set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")
