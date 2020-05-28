include(EkatCreateUnitTest)

# Test kokkos utils
EkatCreateUnitTest(kokkos_utils kokkos_utils_tests.cpp
  LIBS ekat
  THREADS 1 ${EKAT_TEST_MAX_THREADS} ${EKAT_TEST_THREAD_INC}
  INCLUDE_DIRS "${EKAT_INCLUDE_DIRS}"
)

# Test workspace manager
EkatCreateUnitTest(wsm workspace_tests.cpp
  LIBS ekat
  THREADS 1 ${EKAT_TEST_MAX_THREADS} ${EKAT_TEST_THREAD_INC}
  INCLUDE_DIRS "${EKAT_INCLUDE_DIRS}"
)
