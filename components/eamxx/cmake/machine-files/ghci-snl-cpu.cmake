# Common settings for our ghci images
include(${CMAKE_CURRENT_LIST_DIR}/ghci-snl.cmake)

# Set SCREAM_MACHINE
set(SCREAM_MACHINE ghci-snl-cpu CACHE STRING "")

option (EAMXX_ENABLE_PYTHON "Whether to enable python interface from eamxx" ON)
set (Python_EXECUTABLE "/usr/bin/python3" CACHE STRING "")
