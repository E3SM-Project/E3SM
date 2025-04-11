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

if (compile_threaded)
  string(APPEND CMAKE_C_FLAGS " -fopenmp")
  string(APPEND CMAKE_Fortran_FLAGS " -fopenmp")
  string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fopenmp")
endif()

string(APPEND CMAKE_Fortran_FLAGS " -hipa0 -hzero -f free")

string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/opt/gcc/12.2.0/snos/lib64")

# Crusher: this resolves a crash in mct in docn init
string(APPEND CMAKE_C_FLAGS_RELEASE " -O2 -hnoacc -hfp0 -hipa0")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2 -hnoacc -hfp0 -hipa0")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2 ")

string(APPEND CPPDEFS " -DCPRCRAY")

if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()
set(PIO_FILESYSTEM_HINTS "lustre")
