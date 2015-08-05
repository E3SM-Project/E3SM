cmake_minimum_required (VERSION 3.0.0)
project (PIO C)

option (ENABLE_FORTRAN "Enable the Fortran library builds" OFF)
option (ENABLE_TIMING  "Enable GPTL timing in libraries"   OFF)
option (ENABLE_TESTS   "Enable the testing builds"         OFF)
     
add_subdirectory (src)
