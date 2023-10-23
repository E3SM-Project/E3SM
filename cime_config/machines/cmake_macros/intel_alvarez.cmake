string(APPEND CONFIG_ARGS " --host=cray")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
set(CXX_LINKER "FORTRAN")

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
set(CXXFLAGS " ") # hardcode it here to blank, then try to do same things as in intel.cmake
if (compile_threaded)
  string(APPEND CXXFLAGS " -qopenmp")
endif()
if (DEBUG)
  string(APPEND CXXFLAGS " -O0 -g")
endif()
if (NOT DEBUG)
  string(APPEND CXXFLAGS " -O2")
endif()
string(APPEND CXXFLAGS " -fp-model=precise") # and manually add precise
#message(STATUS "ndk CXXFLAGS=${CXXFLAGS}")

string(APPEND FFLAGS " -fp-model=consistent -fimf-use-svml")
if (NOT DEBUG)
   #  string(APPEND FFLAGS " -qno-opt-dynamic-align")
   string(APPEND FFLAGS " -g -traceback")
   string(APPEND CXXFLAGS " -g -traceback")
endif()
string(APPEND FFLAGS " -DHAVE_ERF_INTRINSICS")
string(APPEND CXXFLAGS " -fp-model=consistent")
