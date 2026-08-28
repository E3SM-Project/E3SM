# Get common intel flags
include (${CMAKE_CURRENT_SOURCE_DIR}/intel.cmake)

string(APPEND CMAKE_EXE_LINKER_FLAGS " -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")

string(APPEND CMAKE_CXX_FLAGS_RELEASE " --offload-compress")
