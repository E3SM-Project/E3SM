# Get common intel flags
include (${CMAKE_CURRENT_LIST_DIR}/intel.cmake)

string(APPEND CMAKE_EXE_LINKER_FLAGS " -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")

