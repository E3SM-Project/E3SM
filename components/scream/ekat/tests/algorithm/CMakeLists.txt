include(EkatCreateUnitTest)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/ekat_config.f.in
                ${CMAKE_CURRENT_BINARY_DIR}/ekat_config.f)

# Test lin interp processes
set (LIN_INTERP_SRCS
  lin_interp_test.cpp
  lin_interp_ref.f90
)
EkatCreateUnitTest(lin_interp "${LIN_INTERP_SRCS}"
  LIBS ekat
  THREADS 1 ${EKAT_TEST_MAX_THREADS} ${EKAT_TEST_THREAD_INC}
  INCLUDE_DIRS "${EKAT_INCLUDE_DIRS}")

# Test tridiag solvers
set (TRIDIAG_SRCS
  tridiag_tests.cpp
  tridiag_bfb.cpp
  tridiag_tests_correctness.cpp
  tridiag_tests_performance.cpp
  tridiag_bfb_mod.F90
)
EkatCreateUnitTest(tridiag "${TRIDIAG_SRCS}"
  LIBS ekat
  THREADS ${EKAT_TEST_MAX_THREADS}
  INCLUDE_DIRS "${EKAT_INCLUDE_DIRS}"
  OPTIONAL EXCLUDE_MAIN_CPP)
