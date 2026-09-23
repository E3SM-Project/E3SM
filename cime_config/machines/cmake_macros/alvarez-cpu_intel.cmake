string(APPEND CONFIG_ARGS " --host=cray")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()

set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")

# Check for Intel LLVM (ifx) version 2025 or newer
if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
        #message(STATUS "ndk: Applying Intel 2025.3 sanitization workaround")
        string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check nouninit") # Applying Intel 2025.3 sanitization workaround
    endif()
endif()

#  string(APPEND FFLAGS " -qno-opt-dynamic-align")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -g -traceback")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -g -traceback")

string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fpe0 -init=snan,arrays")
if (COMP_NAME STREQUAL cice)
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -init=nosnan,arrays")
endif()

