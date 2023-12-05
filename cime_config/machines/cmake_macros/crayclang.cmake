if (compile_threaded)
  string(APPEND CMAKE_C_FLAGS   " -fopenmp")
  string(APPEND CMAKE_Fortran_FLAGS   " -fopenmp")
  string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
  string(APPEND CMAKE_EXE_LINKER_FLAGS  " -fopenmp")
endif()
string(APPEND CMAKE_C_FLAGS_DEBUG   " -O0 -g")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g")
string(APPEND CPPDEFS_DEBUG  " -DYAKL_DEBUG")
string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRCRAY")
# -em (default) generates MODULENAME.mod files
string(APPEND CMAKE_Fortran_FLAGS " -f free -em")
if (NOT compile_threaded)
	# -M1077 flag used to suppress message about OpenMP directives
	# that are ignored for non-threaded builds. (-h omp inactive)
	# Details: `explain ftn-1077`
  string(APPEND CMAKE_Fortran_FLAGS " -M1077")
endif()
set(HAS_F2008_CONTIGUOUS "TRUE")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,--allow-multiple-definition")

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

