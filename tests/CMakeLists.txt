cmake_minimum_required (VERSION 3.0.0)
project (PIOTests C Fortran)
enable_testing()
include (CTest)

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================

add_subdirectory (unit)
add_subdirectory (general)
if (PIO_ENABLE_TIMING)
    add_subdirectory (performance)
endif ()
