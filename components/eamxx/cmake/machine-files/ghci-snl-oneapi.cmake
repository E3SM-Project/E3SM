# Common settings for our ghci images
include(${CMAKE_CURRENT_LIST_DIR}/ghci-snl.cmake)

# Set SCREAM_MACHINE
set(SCREAM_MACHINE ghci-snl-oneapi CACHE STRING "")

# Link to MKL
set (CMAKE_EXE_LINKER_FLAGS "-qmkl" CACHE STRING "")

option (EAMXX_ENABLE_PYTHON "Whether to enable python interface from eamxx" ON)
