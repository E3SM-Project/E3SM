string(APPEND CONFIG_ARGS " --host=cray")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()

set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "icx")
set(SCXX "icpx")
set(SFC "ifx")

# Bit of a hack here. For whatever reason, the intel version on pm-cpu (both intel and intel-oneapi, and both icpc/icpx)
# does not seem to have the -fp-model=source flag (docs still show it).  And I was unable to find a reliable way of testing
# on the compiler ID or version, so for now, simply manually adjust the CXXFLAG setting for pm-cpu/intel
# Try to manually remove -fp-model=source (and replace with -fp-model=precise) from CXXFLAGS
#message(STATUS "ndk CXXFLAGS=${CXXFLAGS}")
set(CMAKE_CXX_FLAGS " ") # hardcode it here to blank, then try to do same things as in intel.cmake
if (compile_threaded)
  string(APPEND CMAKE_CXX_FLAGS " -qopenmp")
endif()
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CXX_FLAGS " -fp-model=precise") # and manually add precise
#message(STATUS "ndk CXXFLAGS=${CXXFLAGS}")

string(APPEND CMAKE_Fortran_FLAGS " -fp-model=consistent -fimf-use-svml")
   #  string(APPEND FFLAGS " -qno-opt-dynamic-align")
 string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -g -traceback")
 string(APPEND CMAKE_CXX_FLAGS_RELEASE " -g -traceback")
string(APPEND CMAKE_Fortran_FLAGS " -DHAVE_ERF_INTRINSICS")
string(APPEND CMAKE_CXX_FLAGS " -fp-model=consistent")
