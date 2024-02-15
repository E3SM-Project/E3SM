# Ran into Scorpio build issues using base and child compiler macro files.
# Especially specified compilers were not picked up properly.
# Had to replicate section below to get things working as expected.
# Commenting out until we figure out how to properly handle this. 

# if (compile_threaded)
#   string(APPEND FFLAGS   " -fopenmp")
#   string(APPEND CFLAGS   " -fopenmp")
#   string(APPEND CXXFLAGS " -fopenmp")
#   string(APPEND LDFLAGS  " -fopenmp")
# endif()
# if (DEBUG)
#   string(APPEND CFLAGS   " -O0 -g")
#   string(APPEND FFLAGS   " -O0 -g")
#   string(APPEND CXXFLAGS " -O0 -g")
#   string(APPEND CPPDEFS " -DYAKL_DEBUG")
# endif()
# string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRCRAY")
# string(APPEND FFLAGS " -f free  -em")
# if (NOT compile_threaded)
#   # -M1077 flag used to suppress message about OpenMP directives
#   # that are ignored for non-threaded builds. (-h omp inactive)
#   # Details: `explain ftn-1077`
#   string(APPEND FFLAGS " -M1077")
# endif()
# string(APPEND FFLAGS_NOOPT " -O0")
# set(HAS_F2008_CONTIGUOUS "TRUE")
# string(APPEND LDFLAGS " -Wl,--allow-multiple-definition")
# set(MPICC "cc")
# set(MPICXX "hipcc")
# set(MPIFC "ftn")
# set(SCC "cc")
# set(SCXX "hipcc")
# set(SFC "ftn")

#if (NOT DEBUG)
#  string(APPEND CFLAGS   " -O2")
#  string(APPEND CXXFLAGS " -O2")
#  string(APPEND FFLAGS   " -O2")
#endif()


