string(APPEND CONFIG_ARGS " --host=cray")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -g")
if (compile_threaded)
  string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_OPENMP=Off") # work-around for nvidia as kokkos is not passing "-mp" for threaded build
endif()

# For MCT coupler sources only, try resetting CMAKE_Fortran_FLAGS_DEBUG to avoid setting FPE invalid exceptions
# https://github.com/E3SM-Project/E3SM/issues/7049
if (COMP_NAME STREQUAL cpl)
  set(CMAKE_Fortran_FLAGS_DEBUG " ") # set to blank and rebuild
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -i4 -Mstack_arrays  -Mextend -byteswapio -Mflushz -Kieee -DHAVE_IEEE_ARITHMETIC -Mallocatable=03 -DNO_R16 -traceback")
  #original string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -Ktrap=fp -Mbounds -Kieee")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -Mbounds -Kieee")
  if (compile_threaded)
    string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -mp")
  endif()
endif()

set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")
