include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE)
