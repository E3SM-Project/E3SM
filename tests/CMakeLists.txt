cmake_minimum_required (VERSION 2.8.11)
project (PIOTests C Fortran)

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================

add_subdirectory (unit)
add_subdirectory (general)
if (PIO_ENABLE_TIMING)
    add_subdirectory (performance)
else ()
    message (STATUS "Cannot build performance test without gptl timing library")
endif ()
