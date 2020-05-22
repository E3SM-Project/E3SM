include(EkatCreateUnitTest)

# Test kokkos utils
EkatCreateUnitTest(kokkos_utils kokkos_utils_tests.cpp
  LIBS ekat
  THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
  INCLUDE_DIRS "${EKAT_INCLUDE_DIRS}"
)

# Test workspace manager
EkatCreateUnitTest(wsm workspace_tests.cpp
  LIBS ekat
  THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
  INCLUDE_DIRS "${EKAT_INCLUDE_DIRS}"
)
