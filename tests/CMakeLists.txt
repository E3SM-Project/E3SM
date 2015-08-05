cmake_minimum_required (VERSION 3.0.0)
project (PIOTests C Fortran)

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================

add_subdirectory (general)
#add_subdirectory (unit)
#add_subdirectory (performance)

#==============================================================================
#  DEFINE THE TESTS
#==============================================================================

enable_testing()

#add_test(...)