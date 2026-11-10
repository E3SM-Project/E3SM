include (${CMAKE_CURRENT_LIST_DIR}/gnu.cmake)

string(APPEND CMAKE_Fortran_FLAGS_DEBUG   " -fcheck=bounds")
