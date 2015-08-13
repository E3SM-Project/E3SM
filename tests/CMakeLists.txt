cmake_minimum_required (VERSION 3.0.0)
project (PIOTests C Fortran)
enable_testing()
include (CTest)
include (CheckPlatform)

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================
check_platform ()

add_subdirectory (unit)
add_subdirectory (general)
if (PIO_ENABLE_TIMING)
    add_subdirectory (performance)
else ()
    message (STATUS "Cannot build performance test without gptl timing library")
endif ()
