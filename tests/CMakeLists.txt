cmake_minimum_required (VERSION 2.8.12)
project (PIOTests C Fortran)

#==============================================================================
#  TESTING TARGET
#==============================================================================

# Custom "check" target (builds and runs ctest)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

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
