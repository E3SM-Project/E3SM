include(ScreamUtils)

# Tests to ensure that fails/errors are indeed caught
# NOTE: we don't need any libs for this test, but scream's CreateUnitTest
# has libs as a required arg. So use the raw EkatCreateUnitTest

# This serves as a base case for the following fail checks, since it verifies that
#  - a REQUIRE clause that fails makes the test fail
#  - Our Create unit test logic does work for  catching failures
CreateUnitTest (fail "fail.cpp"
  WILL_FAIL LABELS "fail")

if (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK)
  # Ensure that Kokkos OOB are caught
  CreateUnitTest (kokkos_fail "kokkos_fail.cpp"
    WILL_FAIL LABELS "fail")
endif()

if (EKAT_ENABLE_VALGRIND)
  # Ensure that valgrind errors are caught
  EkatCreateUnitTest (valg_fail "valg_fail.cpp"
    WILL_FAIL LABELS "fail")
endif()

# Ensure that FPE *do* throw when we expect them to
CreateUnitTestExec (scream_fpe_check "fpe_check.cpp")
if (SCREAM_FPE)
  CreateUnitTestFromExec (scream_fpe_check scream_fpe_check
    WILL_FAIL LABELS "check")
else()
  CreateUnitTestFromExec (scream_fpe_check scream_fpe_check
    LABELS "check")
endif()
