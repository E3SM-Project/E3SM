include(ScreamUtils)

# Thin C bindings to individual Fortran COSP routines, used as the BFB reference for
# the C++/Kokkos port. The cosp library writes its .mod files to its own binary dir.
add_library(cosp_f90_shims cosp_f90_shims.F90)
target_link_libraries(cosp_f90_shims PUBLIC cosp)
target_include_directories(cosp_f90_shims PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/..)

# BFB tests of the C++ leaf routines against the Fortran
CreateUnitTest(cosp_leaf_tests
  SOURCES cosp_leaf_tests.cpp
  LIBS eamxx_cosp cosp_f90_shims
  LABELS "cosp;physics")
