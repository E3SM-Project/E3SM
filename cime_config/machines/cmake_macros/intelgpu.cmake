# Get common intel flags
include (${CMAKE_CURRENT_LIST_DIR}/intel.cmake)

# 'just' -g may lead to linker internal errors and/or huge builds out of quotas
string(APPEND CMAKE_C_FLAGS_DEBUG   " -fno-system-debug")
string(APPEND CMAKE_CXX_FLAGS_DEBUG   " -fno-system-debug")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " --offload-compress")

