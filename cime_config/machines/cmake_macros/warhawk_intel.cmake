string(APPEND CMAKE_Fortran_FLAGS " -fimf-use-svml")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -qno-opt-dynamic-align")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lpthread")
