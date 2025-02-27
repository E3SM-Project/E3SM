string(APPEND CONFIG_ARGS " --host=cray")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CMAKE_C_FLAGS_RELEASE " -O2 -g")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2 -g")
set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "gcc")
set(SCXX "g++")
set(SFC "gfortran")

#string(APPEND CMAKE_EXE_LINKER_FLAGS " -static-libstdc++") # was causing link error after Feb 18/19 maintenance

# https://github.com/E3SM-Project/E3SM/issues/7049
if (DEBUG)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " /opt/cray/pe/lib64/libhdf5_hl_parallel_gnu_91.so.200")
endif()
