include(ScreamUtils)

# Thin C bindings to individual Fortran COSP routines, used as the BFB reference for
# the C++/Kokkos port. The cosp library writes its .mod files to its own binary dir.
add_library(cosp_f90_shims cosp_f90_shims.F90)
target_link_libraries(cosp_f90_shims PUBLIC cosp)
target_include_directories(cosp_f90_shims PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/..)

# Exact agreement with the Fortran is only required with the strict floating point model,
# since otherwise the C++ and Fortran compilers may contract and vectorize differently
set(COSP_TEST_DEFS)
string(TOLOWER "${SCREAM_FPMODEL}" COSP_FPMODEL)
if (COSP_FPMODEL STREQUAL "strict")
  list(APPEND COSP_TEST_DEFS COSP_TEST_EXACT_BFB)
endif()

# BFB tests of the C++ leaf routines against the Fortran
CreateUnitTest(cosp_leaf_tests
  SOURCES cosp_leaf_tests.cpp cosp_misr_tests.cpp cosp_isccp_tests.cpp
          cosp_modis_tests.cpp cosp_scops_tests.cpp cosp_main_tests.cpp
  LIBS eamxx_cosp cosp_f90_shims
  COMPILER_CXX_DEFS ${COSP_TEST_DEFS}
  LABELS "cosp;physics")

# The whole COSP call with 10 subcolumns, in its own executable since the Fortran COSP can
# only be initialized once per process
CreateUnitTest(cosp_main_nsubcol10_tests
  SOURCES cosp_main_tests.cpp
  LIBS eamxx_cosp cosp_f90_shims
  COMPILER_CXX_DEFS COSP_TEST_NSUBCOL=10 ${COSP_TEST_DEFS}
  LABELS "cosp;physics")
