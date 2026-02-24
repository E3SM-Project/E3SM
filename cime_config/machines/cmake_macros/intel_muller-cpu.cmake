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


#message(STATUS "ndk CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")
#message(STATUS "ndk CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
#message(STATUS "ndk CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
#message(STATUS "ndk CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
#message(STATUS "ndk CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
#message(STATUS "ndk CMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}")
#-- ndk CMAKE_Fortran_COMPILER_ID=IntelLLVM
#-- ndk CMAKE_Fortran_COMPILER=/opt/cray/pe/craype/2.7.35/bin/ftn
#-- ndk CMAKE_Fortran_COMPILER_VERSION=2025.3.0
#-- ndk CMAKE_CXX_COMPILER_ID=IntelLLVM
#-- ndk CMAKE_CXX_COMPILER=/opt/cray/pe/craype/2.7.35/bin/CC
#-- ndk CMAKE_CXX_COMPILER_VERSION=2025.3.1

# Check for Intel LLVM (ifx) version 2025 or newer
if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
        #message(STATUS "ndk: Applying Intel 2025.3 sanitization workaround")
        string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -init=none -check nouninit") # Applying Intel 2025.3 sanitization workaround (to revert -init=snan,arrays)
        #string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -init=none")
        #string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check nouninit") # may be all that is needed
        #string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fno-sanitize=all") #ndk or maybe just -fno-sanitize=memory ?
        #string(APPEND CMAKE_C_FLAGS_DEBUG " -fno-sanitize=all") #ndk
        #string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fno-sanitize=all") #ndk -- needed?
    endif()
endif()


string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CXX_FLAGS " -fp-model=precise") # and manually add precise

string(APPEND CMAKE_Fortran_FLAGS " -fp-model=consistent -fimf-use-svml")
#  string(APPEND FFLAGS " -qno-opt-dynamic-align")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -g -traceback")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -g -traceback")
string(APPEND CMAKE_Fortran_FLAGS " -DHAVE_ERF_INTRINSICS")
string(APPEND CMAKE_CXX_FLAGS " -fp-model=consistent")
