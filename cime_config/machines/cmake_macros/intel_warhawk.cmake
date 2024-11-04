string(APPEND CMAKE_Fortran_FLAGS " -fp-model consistent -fimf-use-svml")
string(APPEND CMAKE_CXX_FLAGS " -fp-model consistent")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -qno-opt-dynamic-align")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lpthread")
