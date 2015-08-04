cmake_minimum_required (VERSION 3.0.0)
project (PIO C)

set (ENABLE_FORTRAN ON
     CACHE BOOL "Enable the Fortran library builds")

set (ENABLE_TIMING ON
     CACHE BOOL "Enable the GPTL timing library build")

set (ENABLE_TESTS OFF
     CACHE BOOL "Enable the testing builds")
     
add_subdirectory (src)

if (ENABLE_TESTS)
    add_subdirectory (tests)
endif ()
