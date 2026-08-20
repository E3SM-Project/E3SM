include (${CMAKE_CURRENT_SOURCE_DIR}/intel.cmake)

string(APPEND CMAKE_CXX_FLAGS_RELEASE " --offload-compress")
